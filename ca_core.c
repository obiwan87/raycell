#include "ca_core.h"

#include <errno.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Initialize a grid container to a known state.
 * - dim: logical dimension (1 or 2 for current usage).
 * - cellCount: active cell count (<= CELLS).
 *
 * The cell memory is zeroed (so all cells become DEAD if DEAD == 0).
 */
void grid_init_empty(Grid *grid, int dim, int cellCount) {
    grid->dim = dim;
    grid->size = cellCount;
    memset(grid->cells, 0, sizeof(grid->cells));
}

/* Randomize the grid by setting each cell to ALIVE with probability percentAlive%.
 * - percentAlive: 0..100.
 *
 * NOTE: Uses rand(); caller should seed srand().
 */
void grid_randomize(Grid *grid, int percentAlive) {
    for (int i = 0; i < grid->size; i++) {
        grid->cells[i].state = (rand() % 100 < percentAlive) ? CELL_STATE_ALIVE : CELL_STATE_DEAD;
    }
}

/* Clone grid metadata and cell states from srcGrid to dstGrid. */
void grid_clone(const Grid *srcGrid, Grid *dstGrid) {
    memcpy(dstGrid->cells, srcGrid->cells, sizeof(Cell) * (size_t)srcGrid->size);
    dstGrid->size = srcGrid->size;
    dstGrid->dim  = srcGrid->dim;
}

/* Apply Wolfram [L,C,R] rule to produce next cells from current cells.
 * - currentGrid: current generation (1D).
 * - ruleCode: Wolfram rule code (0..255).
 * - outCells: caller-provided output buffer (length currentGrid->size).
 *
 * Boundary condition:
 * - The first and last cell are forced to DEAD (fixed zero boundary).
 */
void wolfram_apply_rule_1d(const Grid *currentGrid, WolframCode ruleCode, Cell *outCells) {
    outCells[0].state = CELL_STATE_DEAD;
    outCells[currentGrid->size - 1].state = CELL_STATE_DEAD;

    for (int i = 1; i < currentGrid->size - 1; i++) {
        const int l = currentGrid->cells[i - 1].state;
        const int c = currentGrid->cells[i].state;
        const int r = currentGrid->cells[i + 1].state;

        const int pattern = (l << 2) | (c << 1) | (r << 0);
        outCells[i].state = (ruleCode >> pattern) & 1;
    }
}

/* Initialize a 1D Wolfram run with an empty grid at step 0. */
void wolfram_run_init(WolframRun *run, WolframCode ruleCode) {
    run->current_step = 0;
    run->code = ruleCode;
    grid_init_empty(&run->grids[0], 1, CELLS);
}

/* Advance the Wolfram run by one step (if there is room in grids[]). */
void wolfram_run_step(WolframRun *run) {
    if (run->current_step + 1 >= MAX_STEPS) return;

    const Grid *currentGrid = &run->grids[run->current_step];
    Grid *nextGrid = &run->grids[run->current_step + 1];

    nextGrid->dim  = currentGrid->dim;
    nextGrid->size = currentGrid->size;

    wolfram_apply_rule_1d(currentGrid, run->code, nextGrid->cells);
    run->current_step++;
}

/* Configure Von Neumann neighborhood (4-neighbors):
 * left, right, up, down.
 */
void van_neumann_neighborhood(Neighborhood* neighborhood) {
    neighborhood->size = 8;
    neighborhood->offsets[0] = -1; neighborhood->offsets[1] =  0;
    neighborhood->offsets[2] =  1; neighborhood->offsets[3] =  0;
    neighborhood->offsets[4] =  0; neighborhood->offsets[5] = -1;
    neighborhood->offsets[6] =  0; neighborhood->offsets[7] =  1;
}

/* Configure Moore neighborhood (8-neighbors):
 * the 8 surrounding cells.
 *
 * NOTE:
 * There is a duplicated (-1, +1) and missing (-1, -1) in the original listing.
 * I am preserving your code as-is, only documenting. Fix if unintended.
 */
void moore_neighborhood(Neighborhood* neighborhood) {
    neighborhood->size = 16;
    neighborhood->offsets[0]  = -1; neighborhood->offsets[1]  =  1;
    neighborhood->offsets[2]  = -1; neighborhood->offsets[3]  =  0;
    neighborhood->offsets[4]  = -1; neighborhood->offsets[5]  =  1;
    neighborhood->offsets[6]  =  0; neighborhood->offsets[7]  = -1;
    neighborhood->offsets[8]  =  0; neighborhood->offsets[9]  =  1;
    neighborhood->offsets[10] =  1; neighborhood->offsets[11] = -1;
    neighborhood->offsets[12] =  1; neighborhood->offsets[13] =  0;
    neighborhood->offsets[14] =  1; neighborhood->offsets[15] =  1;
}

/* Count how many neighbor cells are ALIVE in a provided neighbor list.
 * - neighborCells: array of pointers (some entries may be NULL depending on topology).
 * - neighborCount: number of pointers in neighborCells.
 *
 * Returns:
 * - number of neighbors with state == ALIVE.
 */
int ca_count_alive(Cell** neighborCells, size_t neighborCount) {
    int alive = 0;
    for (size_t i = 0; i < neighborCount; i++) {
        if (neighborCells[i] && neighborCells[i]->state == CELL_STATE_ALIVE) {
            alive++;
        }
    }
    return alive;
}

/* Fill out_neighbors with pointers to the neighbors of the cell at `cellIndex`.
 *
 * Topology:
 * - Torus (wrap-around) via modulo arithmetic on x/y coordinates.
 *
 * Inputs:
 * - automaton: provides current grid and neighborhood offsets.
 * - cellIndex: index in row-major order.
 * Output:
 * - out_neighbors: must have capacity neighborhood.size/2.
 */
void ca_get_neighbors(CellularAutomaton* automaton, int cellIndex, Cell** out_neighbors) {
    const int side = (int)sqrt((double)automaton->currentGenBuf.size);

    /* Convert linear index to (x,y). */
    const int x = cellIndex % side;
    const int y = cellIndex / side;

    /* For each offset pair, compute wrapped neighbor coordinates. */
    int out_i = 0;
    for (int i = 0; i + 1 < automaton->neighborhood.size; i += 2) {
        const int dx = automaton->neighborhood.offsets[i];
        const int dy = automaton->neighborhood.offsets[i + 1];

        int nx = (x + dx) % side;
        int ny = (y + dy) % side;
        if (nx < 0) nx += side;
        if (ny < 0) ny += side;

        const int n_idx = ny * side + nx;
        out_neighbors[out_i++] = &automaton->currentGenBuf.cells[n_idx];
    }
}

/* Apply a generic B/S rule:
 * - b: digits for birth counts (when current is DEAD)
 * - s: digits for survival counts (when current is ALIVE)
 *
 * Example:
 * - Game of Life is B3/S23
 *
 * Returns:
 * - next state (ALIVE or DEAD).
 *
 * NOTE:
 * This parser assumes only single-digit neighbor counts, which is correct for
 * Moore/VonNeumann with <= 8 neighbors. If you add larger neighborhoods, you’ll
 * need a different encoding (or multi-digit parsing).
 */
int ca_apply_bs_rule(const char* b, const char* s, int currentState, int aliveNeighbors) {
    if (currentState == CELL_STATE_DEAD) {
        const int b_len = (int)strlen(b);
        bool birth = false;
        for (int i = 0; i < b_len && !birth; i++) {
            const int n = b[i] - '0';
            birth = (aliveNeighbors == n);
        }
        return birth ? CELL_STATE_ALIVE : CELL_STATE_DEAD;
    }

    if (currentState == CELL_STATE_ALIVE) {
        const int s_len = (int)strlen(s);
        bool survival = false;
        for (int i = 0; i < s_len && !survival; i++) {
            const int n = s[i] - '0';
            survival = (aliveNeighbors == n);
        }
        return survival ? CELL_STATE_ALIVE : CELL_STATE_DEAD;
    }

    /* Fallback for unexpected states. */
    return CELL_STATE_DEAD;
}

/* Transition function: Conway's Game of Life (B3/S23).
 * - env is unused.
 */
int ca_transition_game_of_life(const void* env, const Cell* currentCell, Cell** neighbors, size_t neighborCount) {
    (void) env;
    return ca_apply_bs_rule("3", "23", currentCell->state, ca_count_alive(neighbors, neighborCount));
}

/* Transition function: "Replicator" rule (B1357/S1357).
 * - env is unused.
 */
int ca_transition_replicator(const void* env, const Cell* currentCell, Cell** neighbors, size_t neighborCount) {
    (void) env;
    return ca_apply_bs_rule("1357", "1357", currentCell->state, ca_count_alive(neighbors, neighborCount));
}

/* Environment payload for a generic B/S rule transition. */
typedef struct {
    char* b; /* Birth digits string, e.g. "3" or "1357". */
    char* s; /* Survival digits string, e.g. "23" or "1357". */
} BSRule;

/* Transition function: generic B/S rule using BSRule env payload. */
static int ca_transition_bs_rule(const void* env, const Cell* currentCell, Cell** neighbors, size_t neighborCount) {
    const BSRule* bsRule = (const BSRule*) env;
    return ca_apply_bs_rule(bsRule->b, bsRule->s, currentCell->state, ca_count_alive(neighbors, neighborCount));
}

/* Duplicate a C-string into a freshly allocated buffer.
 * Returns NULL on error and sets errno to EINVAL if src == NULL.
 */
static char* dup_cstr_or_null(const char *src) {
    if (!src) {
        errno = EINVAL;
        return NULL;
    }
    const size_t n = strlen(src);
    char *dst = (char*)malloc(n + 1);
    if (!dst) return NULL;
    memcpy(dst, src, n + 1); /* includes '\0' */
    return dst;
}

/* Allocate a Transition configured with a heap-allocated BSRule env payload.
 * The returned Transition must be freed via transition_free_bs_rule().
 *
 * Inputs:
 * - b: Birth digits string (e.g. "3")
 * - s: Survival digits string (e.g. "23")
 *
 * Returns:
 * - Transition* on success, NULL on error (errno may be set).
 */
Transition* transition_alloc_bs_rule(const char *b, const char *s) {
    if (b == NULL || s == NULL) {
        errno = EINVAL;
        return NULL;
    }

    Transition *transition = (Transition*)malloc(sizeof *transition);
    if (!transition) return NULL;

    BSRule *bsRule = (BSRule*)malloc(sizeof *bsRule);
    if (!bsRule) {
        free(transition);
        return NULL;
    }

    bsRule->b = dup_cstr_or_null(b);
    if (!bsRule->b) {
        free(bsRule);
        free(transition);
        return NULL;
    }

    bsRule->s = dup_cstr_or_null(s);
    if (!bsRule->s) {
        free(bsRule->b);
        free(bsRule);
        free(transition);
        return NULL;
    }

    transition->env = bsRule;
    transition->fn = &ca_transition_bs_rule;

    return transition;
}

/* Free a Transition previously allocated with transition_alloc_bs_rule(). */
void transition_free_bs_rule(Transition* t) {
    if (t == NULL) return;

    BSRule* bsRule = (BSRule*)t->env;
    free(bsRule->b);
    free(bsRule->s);
    free(bsRule);
    free(t);
}

/* Grid initializer implementation:
 * Randomly marks percentAlive% of the cells as ALIVE.
 *
 * env is expected to be int* (percentAlive).
 */
static void grid_init_random_percent_alive(const void* env, Grid* grid) {
    const int percentAlive = *(const int*)env;
    grid_randomize(grid, percentAlive);
}

/* Allocate a GridInitializer that seeds the grid with percentAlive% alive cells.
 * The returned initializer must be freed via grid_initializer_free_percent_alive().
 */
GridInitializer* grid_initializer_alloc_percent_alive(int percentAlive) {
    GridInitializer* gridInitializer = (GridInitializer*)malloc(sizeof(GridInitializer));
    if(!gridInitializer) return NULL;

    int* pPercentAlive = (int*)malloc(sizeof(int));
    if(!pPercentAlive) {
        free(gridInitializer);
        return NULL;
    }

    *pPercentAlive = percentAlive;

    gridInitializer->env = pPercentAlive;
    gridInitializer->fn = &grid_init_random_percent_alive;

    return gridInitializer;
}

/* Free a GridInitializer previously allocated with grid_initializer_alloc_percent_alive(). */
void grid_initializer_free_percent_alive(GridInitializer* gridInitializer) {
    if (!gridInitializer) return;
    free(gridInitializer->env);
    free(gridInitializer);
}

/* Grid initializer implementation:
 * Clears the grid and sets exactly one random cell to ALIVE.
 *
 * env is unused.
 */
static void grid_init_random_single_alive_cell(const void* env, Grid* grid) {
    (void) env;

    const size_t randomCell = (size_t)rand() % (size_t)grid->size;
    printf("Random Cell: %lu\n", (unsigned long)randomCell);

    for(int i = 0; i < grid->size; i++) {
        grid->cells[i].state = CELL_STATE_DEAD;
    }

    grid->cells[randomCell].state = CELL_STATE_ALIVE;
}

/* Allocate a GridInitializer that seeds the grid with a single alive cell.
 * The returned initializer must be freed via grid_initializer_free_single_cell_alive().
 *
 * NOTE: This initializer has no env payload.
 */
GridInitializer* grid_initializer_alloc_single_cell_alive(void) {
    GridInitializer* gridInitializer = (GridInitializer*)malloc(sizeof(GridInitializer));
    if(!gridInitializer) return NULL;

    gridInitializer->env = NULL;
    gridInitializer->fn = &grid_init_random_single_alive_cell;

    return gridInitializer;
}

/* Free a GridInitializer previously allocated with grid_initializer_alloc_single_cell_alive(). */
void grid_initializer_free_single_cell_alive(GridInitializer* gridInitializer) {
    free(gridInitializer);
}

/* Advance the 2D cellular automaton by one generation.
 *
 * Behavior:
 * 1) For each cell:
 *    - gather neighbor pointers (according to neighborhood offsets)
 *    - call transition function to compute next state
 *    - write into nextGenBuf
 * 2) Save current generation into ring-buffer history (for stepping backwards)
 * 3) Copy nextGenBuf into currentGenBuf
 * 4) Increment historyIndex (acts as generation counter and ring cursor)
 */
void ca_step_forward(CellularAutomaton* ca) {
    /* Neighbor pointer scratch array.
     * Capacity = MAX_NEIGHBORHOOD/2 because offsets are stored as dx,dy pairs.
     */
    Cell* neighbors[MAX_NEIGHBORHOOD / 2] = {0};

    const size_t neighborCount = (size_t)ca->neighborhood.size / 2;

    for (int i = 0; i < ca->currentGenBuf.size; i++) {
        const Cell* currentCell = &ca->currentGenBuf.cells[i];
        ca_get_neighbors(ca, i, neighbors);

        const int newState = ca->transition.fn(ca->transition.env, currentCell, neighbors, neighborCount);
        ca->nextGenBuf.cells[i].state = newState;
    }

    /* Store current generation into history ring buffer. */
    memcpy(
        ca->history + (ca->historyIndex % HISTORY_SIZE) * ca->currentGenBuf.size,
        ca->currentGenBuf.cells,
        (size_t)ca->currentGenBuf.size * sizeof(Cell)
    );

    /* Commit next generation. */
    memcpy(
        ca->currentGenBuf.cells,
        ca->nextGenBuf.cells,
        (size_t)ca->currentGenBuf.size * sizeof(Cell)
    );
    ca->historyIndex++;
}

/* Step the automaton one generation backwards (if history is available).
 * Restores currentGenBuf from the history ring buffer.
 */
void ca_step_back(CellularAutomaton* automaton) {
    if (automaton->historyIndex <= 0) return;

    const int prevIndex = automaton->historyIndex - 1;
    memcpy(
        automaton->currentGenBuf.cells,
        automaton->history + (prevIndex % HISTORY_SIZE) * automaton->currentGenBuf.size,
        (size_t)automaton->currentGenBuf.size * sizeof(Cell)
    );
    automaton->historyIndex--;
}


CellularAutomaton* ca_alloc_2d(const Neighborhood neighborhood)
{
    CellularAutomaton* ca = malloc(sizeof(CellularAutomaton));
    ca->neighborhood = neighborhood;

    ca->currentGenBuf.dim = 2;
    ca->currentGenBuf.size = CELLS;

    ca->nextGenBuf.dim = 2;
    ca->nextGenBuf.size = CELLS;
    ca->historyIndex = 0;
    return ca;
}

void ca_free_2d(CellularAutomaton* ca)
{
    free(ca);
}