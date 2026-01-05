#include "raylib.h"
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

/* Window dimensions in pixels. */
#define WINDOW_WIDTH  800
#define WINDOW_HEIGHT 800

/* Visual size of a single cell in pixels. */
#define CELL_WIDTH   4
#define CELL_HEIGHT  4

/* Total number of cells in the grid.
 * NOTE: This code assumes a square grid where side = sqrt(CELLS).
 */
#define CELLS        (200*200)

/* Canonical binary states used by current rules/rendering.
 * NOTE: The code is mostly prepared for extension, but many rules assume DEAD/ALIVE.
 */
#define DEAD   0
#define ALIVE  1

/* Maximum number of ints stored in Neighborhood.offsets.
 * Because offsets are stored as dx,dy pairs, this corresponds to MAX_NEIGHBORHOOD/2 neighbors.
 */
#define MAX_NEIGHBORHOOD 16

/* Maximum number of steps stored for the 1D Wolfram run (unused in raycell_loop). */
#define MAX_STEPS        2048

/* Number of past generations stored for stepping backwards in the 2D CA demo. */
#define HISTORY_SIZE     10

/* Background color used for clearing the window. */
#define BACKGROUND ((Color){ 50, 50, 50, 255 })

/* A single cell in the automaton grid.
 * - state: integer state. Currently binary (DEAD/ALIVE), but can be extended for multi-state rules.
 */
typedef struct {
    int state;
} Cell;

/* A generic grid container that holds a flat array of cells.
 * - cells: fixed-capacity storage (CELLS).
 * - dim:  logical dimensionality (e.g. 1 for Wolfram, 2 for 2D CA).
 * - size: number of active cells (<= CELLS). For 2D this is expected to be side*side.
 *
 * NOTE:
 * The cells are stored in row-major order for 2D:
 *   index = y * side + x
 */
typedef struct {
    Cell cells[CELLS];
    int dim;
    int size;
} Grid;

/* Wolfram elementary automaton rule number (0..255). */
typedef int WolframCode;

/* A run (history buffer) for a 1D Wolfram automaton visualization.
 * - current_step: index of the last computed generation in grids[].
 * - code: Wolfram rule code (0..255).
 * - grids: array of consecutive generations (0..current_step are valid).
 *
 * NOTE: This is currently unused in raycell_loop(), but kept for experimentation.
 */
typedef struct {
    int current_step;
    WolframCode code;
    Grid grids[MAX_STEPS];
} WolframRun;

/* Neighborhood definition for 2D cellular automata.
 * - offsets: interleaved dx,dy pairs: [dx0,dy0, dx1,dy1, ...].
 * - size: number of ints in offsets[] (must be even).
 *
 * Example:
 *   size = 8 -> 4 neighbors (Von Neumann)
 *   size = 16 -> 8 neighbors (Moore)
 */
typedef struct {
    int size;
    int offsets[MAX_NEIGHBORHOOD];
} Neighborhood;

/* Transition function signature for computing a cell's next state.
 * - env:             user-provided environment pointer (rule parameters, etc.).
 * - currentCell:     pointer to the current cell.
 * - neighboringCells: array of pointers to neighbor cells (length = neighbors).
 * - neighbors:       number of neighbor pointers.
 *
 * Returns:
 * - next state for the current cell.
 */
typedef int  (*TransitionFn)(const void* env, const Cell* currentCell, Cell** neighboringCells, size_t neighbors);

/* Grid initializer function signature.
 * - env: user-provided environment pointer (initializer parameters).
 * - grid: grid to initialize in-place.
 */
typedef void (*GridInitializerFn)(const void *env, Grid* grid);

/* A grid initializer object.
 * - env: opaque parameter storage for initializer logic.
 * - fn:  initializer implementation.
 */
typedef struct {
    void* env;
    GridInitializerFn fn;
} GridInitializer;

/* A transition rule object.
 * - env: opaque parameter storage for transition logic.
 * - fn:  transition implementation.
 */
typedef struct {
    void* env;
    TransitionFn fn;
} Transition;

/* A single automaton configuration slot, combining:
 * - a Transition rule
 * - a GridInitializer for seeding/resetting the grid
 *
 * NOTE: This is used by the "digit-to-select-config" logic in raycell_loop().
 */
typedef struct {
    Transition* transition;
    GridInitializer* gridInitializer;
} CelluarAutomatonConfig;

/* A 2D cellular automaton simulation state.
 * - currentGenBuf: current generation grid.
 * - nextGenBuf:    scratch buffer for computing next generation.
 * - transition:    active rule used to compute next generation.
 * - neighborhood:  neighbor coordinate offsets used for neighbor lookup.
 * - historyIndex:  number of steps performed (also used as the ring-buffer cursor).
 * - history:       ring-buffer of previous cell states for stepping backwards.
 *
 * NOTE:
 * history stores HISTORY_SIZE copies of the entire grid, in a flat array:
 *   history[ (slot * CELLS) + i ]
 */
typedef struct {
    Grid currentGenBuf;
    Grid nextGenBuf;
    Transition transition;
    Neighborhood neighborhood;
    int historyIndex;
    Cell history[CELLS * HISTORY_SIZE];
} CellularAutomaton;

/* Map a cell state to a render color.
 * NOTE: Currently only handles DEAD/ALIVE. If you add multi-state,
 * update this function (e.g. palette lookup).
 */
static inline Color cell_color(int state) {
    return (state == ALIVE) ? (Color){ 180, 80, 10, 255 }
    : (Color){  50, 50, 50, 255 };
}

/* Initialize a grid container to a known state.
 * - dim: logical dimension (1 or 2 for current usage).
 * - cellCount: active cell count (<= CELLS).
 *
 * The cell memory is zeroed (so all cells become DEAD if DEAD == 0).
 */
static void grid_init_empty(Grid *grid, int dim, int cellCount) {
    grid->dim = dim;
    grid->size = cellCount;
    memset(grid->cells, 0, sizeof(grid->cells));
}

/* Randomize the grid by setting each cell to ALIVE with probability percentAlive%.
 * - percentAlive: 0..100.
 *
 * NOTE: Uses rand(); caller should seed srand().
 */
static void grid_randomize(Grid *grid, int percentAlive) {
    for (int i = 0; i < grid->size; i++) {
        grid->cells[i].state = (rand() % 100 < percentAlive) ? ALIVE : DEAD;
    }
}

/* Clone grid metadata and cell states from srcGrid to dstGrid. */
static void grid_clone(const Grid *srcGrid, Grid *dstGrid) {
    memcpy(dstGrid->cells, srcGrid->cells, sizeof(Cell) * srcGrid->size);
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
static void wolfram_apply_rule_1d(const Grid *currentGrid, WolframCode ruleCode, Cell *outCells) {
    outCells[0].state = DEAD;
    outCells[currentGrid->size - 1].state = DEAD;

    for (int i = 1; i < currentGrid->size - 1; i++) {
        const int l = currentGrid->cells[i - 1].state;
        const int c = currentGrid->cells[i].state;
        const int r = currentGrid->cells[i + 1].state;

        const int pattern = (l << 2) | (c << 1) | (r << 0);
        outCells[i].state = (ruleCode >> pattern) & 1;
    }
}

/* Initialize a 1D Wolfram run with an empty grid at step 0. */
static void wolfram_run_init(WolframRun *run, WolframCode ruleCode) {
    run->current_step = 0;
    run->code = ruleCode;
    grid_init_empty(&run->grids[0], 1, CELLS);
}

/* Advance the Wolfram run by one step (if there is room in grids[]). */
static void wolfram_run_step(WolframRun *run) {
    if (run->current_step + 1 >= MAX_STEPS) return;

    const Grid *currentGrid = &run->grids[run->current_step];
    Grid *nextGrid = &run->grids[run->current_step + 1];

    nextGrid->dim  = currentGrid->dim;
    nextGrid->size = currentGrid->size;

    wolfram_apply_rule_1d(currentGrid, run->code, nextGrid->cells);
    run->current_step++;
}

/* Render/update loop for a 1D Wolfram run visualization.
 * - camera: camera to use for scrolling the history vertically.
 * - intervalSeconds: auto-step interval.
 * - timerSeconds: pointer to an accumulator (maintained by caller).
 *
 * Controls:
 * - SPACE: step once
 * - R: reset to step 0
 *
 * NOTE: This function is not used in raycell_loop(), but remains as a separate demo path.
 */
static void wolfram_run_render(WolframRun* run, Camera2D camera, float intervalSeconds, float* timerSeconds) {
    const float dt = GetFrameTime();

    if (IsKeyDown(KEY_SPACE)) {
        wolfram_run_step(run);
    }

    if (IsKeyPressed(KEY_R)) {
        run->current_step = 0;
    }

    *timerSeconds += dt;
    if (*timerSeconds >= intervalSeconds) {
        *timerSeconds -= intervalSeconds;
        wolfram_run_step(run);
    }

    const int screen_h = GetScreenHeight();
    const int grid_pixel_h = (run->current_step + 1) * CELL_HEIGHT;

    /* Auto-scroll once the grid becomes taller than a threshold. */
    if (grid_pixel_h * 4 > 3 * screen_h) {
        camera.target = (Vector2){ 0.0f, (float)grid_pixel_h - 3.0f * (float)screen_h / 4.0f };
    } else {
        camera.target = (Vector2){ 0.0f, 0.0f };
    }

    BeginDrawing();
    ClearBackground(BACKGROUND);

    BeginMode2D(camera);
    for (int step = 0; step <= run->current_step; step++) {
        const Grid* grid = &run->grids[step];
        for (int cell = 0; cell < grid->size; cell++) {
            DrawRectangle(cell * CELL_WIDTH,
                          step * CELL_HEIGHT,
                          CELL_WIDTH, CELL_HEIGHT,
                          cell_color(grid->cells[cell].state));
        }
    }
    EndMode2D();

    EndDrawing();
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
static int ca_count_alive(Cell** neighborCells, size_t neighborCount) {
    int alive = 0;
    for (size_t i = 0; i < neighborCount; i++) {
        if (neighborCells[i] && neighborCells[i]->state == ALIVE) {
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
static void ca_get_neighbors(CellularAutomaton* automaton, int cellIndex, Cell** out_neighbors) {
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
static int ca_apply_bs_rule(const char* b, const char* s, int currentState, int aliveNeighbors) {
    if (currentState == DEAD) {
        const int b_len = (int)strlen(b);
        bool birth = false;
        for (int i = 0; i < b_len && !birth; i++) {
            const int n = b[i] - '0';
            birth = (aliveNeighbors == n);
        }
        return birth ? ALIVE : DEAD;
    }

    if (currentState == ALIVE) {
        const int s_len = (int)strlen(s);
        bool survival = false;
        for (int i = 0; i < s_len && !survival; i++) {
            const int n = s[i] - '0';
            survival = (aliveNeighbors == n);
        }
        return survival ? ALIVE : DEAD;
    }

    /* Fallback for unexpected states. */
    return DEAD;
}

/* Transition function: Conway's Game of Life (B3/S23).
 * - env is unused.
 */
static int ca_transition_game_of_life(const void* env, const Cell* currentCell, Cell** neighbors, size_t neighborCount) {
    (void) env;
    return ca_apply_bs_rule("3", "23", currentCell->state, ca_count_alive(neighbors, neighborCount));
}

/* Transition function: "Replicator" rule (B1357/S1357).
 * - env is unused.
 */
static int ca_transition_replicator(const void* env, const Cell* currentCell, Cell** neighbors, size_t neighborCount) {
    (void) env;
    return  ca_apply_bs_rule("1357", "1357", currentCell->state, ca_count_alive(neighbors, neighborCount));
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
    const size_t n = strlen(src);
    char *dst = (char*)malloc(n + 1);
    if (!dst) return NULL;
    memcpy(dst, src, n + 1); // includes '\0'
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
    if (t == NULL)
        return;
    BSRule* bsRule = t->env;
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
void grid_init_random_percent_alive(const void* env, Grid* grid) {
    const int percentAlive = *(int*)env;
    grid_randomize(grid, percentAlive);
}

/* Allocate a GridInitializer that seeds the grid with percentAlive% alive cells.
 * The returned initializer must be freed via grid_initializer_free_percent_alive().
 */
GridInitializer* grid_initializer_alloc_percent_alive(int percentAlive) {
    GridInitializer* gridInitializer = (GridInitializer*)malloc(sizeof(GridInitializer));
    if(!gridInitializer)
        return NULL;

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
    free(gridInitializer->env);
    free(gridInitializer);
}

/* Grid initializer implementation:
 * Clears the grid and sets exactly one random cell to ALIVE.
 *
 * env is unused.
 */
void grid_init_random_single_alive_cell(const void* env, Grid* grid) {
    (void) env;

    const size_t randomCell = (size_t)rand() % (size_t)grid->size;
    printf("Random Cell: %lu\n", (unsigned long)randomCell);

    for(int i = 0; i < grid->size; i++) {
        grid->cells[i].state = DEAD;
    }

    grid->cells[randomCell].state = ALIVE;
}

/* Allocate a GridInitializer that seeds the grid with a single alive cell.
 * The returned initializer must be freed via grid_initializer_free_single_cell_alive().
 *
 * NOTE: This initializer has no env payload.
 */
GridInitializer* grid_initializer_alloc_single_cell_alive() {
    GridInitializer* gridInitializer = (GridInitializer*)malloc(sizeof(GridInitializer));
    if(!gridInitializer)
        return NULL;

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
static void ca_step_forward(CellularAutomaton* ca) {
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
    memcpy(ca->currentGenBuf.cells, ca->nextGenBuf.cells, (size_t)ca->currentGenBuf.size * sizeof(Cell));
    ca->historyIndex++;
}

/* Step the automaton one generation backwards (if history is available).
 * Restores currentGenBuf from the history ring buffer.
 */
static void ca_step_back(CellularAutomaton* automaton) {
    if (automaton->historyIndex <= 0) return;

    const int prevIndex = automaton->historyIndex - 1;
    memcpy(
        automaton->currentGenBuf.cells,
        automaton->history + (prevIndex % HISTORY_SIZE) * automaton->currentGenBuf.size,
        (size_t)automaton->currentGenBuf.size * sizeof(Cell)
    );
    automaton->historyIndex--;
}

/* Draw the current 2D automaton grid.
 * Assumptions:
 * - The grid is square.
 * - Cells are stored in a 1D array in row-major order.
 * - size == side*side, where side = sqrt(size).
 */
static void ca_draw_matrix(CellularAutomaton* automaton) {
    const int side = (int)sqrt((double)automaton->currentGenBuf.size);

    for (int row = 0; row < side; row++) {
        for (int col = 0; col < side; col++) {
            const Cell* cell = &automaton->currentGenBuf.cells[row * side + col];
            DrawRectangle(col * CELL_WIDTH,
                          row * CELL_HEIGHT,
                          CELL_WIDTH, CELL_HEIGHT,
                          cell_color(cell->state));
        }
    }
}

/* Global configuration slots:
 * - configs: array of selectable CA configurations.
 * - configs_count: number of initialized entries in configs[].
 * - active_config: currently selected config index.
 *
 * NOTE: Kept global to keep the demo wiring simple.
 */
static CelluarAutomatonConfig configs[10];
static int configs_count = 0;
static int active_config = 0;

/* Main interactive loop for the 2D automaton demo.
 *
 * Controls:
 * - SPACE: step forward once
 * - LEFT:  step backward (within HISTORY_SIZE)
 * - R:     reset current configuration (re-run its initializer)
 * - 0..9:  switch to configuration N (if present) and reset
 */
static void raycell_loop(void) {
    InitWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "raycell");
    SetTargetFPS(60);

    /* Build neighborhood definition used by all configurations. */
    Neighborhood neighborhood;
    moore_neighborhood(&neighborhood);

    /* Config 0: Game of Life (B3/S23), seeded with 40% alive cells. */
    Transition* gofTransition = transition_alloc_bs_rule("3", "23");
    GridInitializer* gofGridInitializer = grid_initializer_alloc_percent_alive(40);
    configs[0].transition = gofTransition;
    configs[0].gridInitializer = gofGridInitializer;
    configs_count++;

    /* Config 1: Replicator (B1357/S1357), seeded with a single alive cell. */
    Transition* replicatorTransition = transition_alloc_bs_rule("1357", "1357");
    GridInitializer* replicatorInitializer = grid_initializer_alloc_single_cell_alive();
    configs[1].transition = replicatorTransition;
    configs[1].gridInitializer = replicatorInitializer;
    configs_count++;

    /* Cellular automaton runtime state. */
    CellularAutomaton* ca = malloc(sizeof(CellularAutomaton));
    ca->neighborhood = neighborhood;

    ca->currentGenBuf.dim = 2;
    ca->currentGenBuf.size = CELLS;

    ca->nextGenBuf.dim = 2;
    ca->nextGenBuf.size = CELLS;
    ca->historyIndex = 0;

    /* Seed initial grid using config 0 initializer. */
    gofGridInitializer->fn(gofGridInitializer->env, &ca->currentGenBuf);

    while (!WindowShouldClose()) {
        /* Detect digit keys (top row and keypad) to switch configs. */
        int digit = -1;
        for (int i = 0; i <= 9; i++) {
            if (IsKeyPressed(KEY_ZERO + i)) { digit = i; break; }
        }
        for (int i = 0; i <= 9 && digit == -1; i++) {
            if (IsKeyPressed(KEY_KP_0 + i)) { digit = i; break; }
        }

        /* Reset logic: either 'R' or selecting a new config. */
        bool shouldReset = IsKeyPressed(KEY_R);
        if(digit >= 0 && digit < configs_count) {
            printf("Digit pressed: %d\n", digit);
            active_config = digit;
            shouldReset = true;
        }

        /* Apply active config to CA. */
        const GridInitializer* activeInitializer = configs[active_config].gridInitializer;
        const Transition* activeTransition = configs[active_config].transition;
        ca->transition = *activeTransition;

        if (shouldReset) {
            activeInitializer->fn(activeInitializer->env, &ca->currentGenBuf);
            ca->historyIndex = 0;
        }

        /* Manual stepping. */
        if (IsKeyPressed(KEY_SPACE)) {
            ca_step_forward(ca);
        }

        if (IsKeyPressed(KEY_LEFT)) {
            ca_step_back(ca);
        }

        /* Draw current generation. */
        BeginDrawing();
        ClearBackground(BACKGROUND);
        ca_draw_matrix(ca);
        EndDrawing();
    }

    CloseWindow();

    /* Cleanup heap-allocated configuration objects. */
    transition_free_bs_rule(gofTransition);
    transition_free_bs_rule(replicatorTransition);

    grid_initializer_free_percent_alive(gofGridInitializer);
    grid_initializer_free_single_cell_alive(replicatorInitializer);

    free(ca);
}

/* Program entry point. */
int main(void) {
    /* Seed PRNG for grid randomization and single-cell seed selection. */
    srand(696969);

    raycell_loop();
    return 0;
}
