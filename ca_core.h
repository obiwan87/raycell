#ifndef CA_CORE_H
#define CA_CORE_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {


#endif

/* Total number of cells in the grid.
* NOTE: This code assumes a square grid where side = sqrt(CELLS).
*/
#define CELLS        (200*200)

/* Canonical binary states used by current rules/rendering.
 * NOTE: The code is mostly prepared for extension, but many rules assume DEAD/ALIVE.
 */
#define CELL_STATE_DEAD   0
#define CELL_STATE_ALIVE  1

/* Maximum number of ints stored in Neighborhood.offsets.
 * Because offsets are stored as dx,dy pairs, this corresponds to MAX_NEIGHBORHOOD/2 neighbors.
 */
#define MAX_NEIGHBORHOOD 16

/* Maximum number of steps stored for the 1D Wolfram run (unused in raycell_loop). */
#define MAX_STEPS        2048

/* Number of past generations stored for stepping backwards in the 2D CA demo. */
#define HISTORY_SIZE     10

/* A single cell in the automaton grid.
 * - state: integer state. Currently binary (DEAD/ALIVE), but can be extended for multi-state rules.
 */
typedef struct
{
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
typedef struct
{
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
typedef struct
{
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
typedef struct
{
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
typedef int (*TransitionFn)(const void* env, const Cell* currentCell, Cell** neighboringCells, size_t neighbors);

/* Grid initializer function signature.
 * - env: user-provided environment pointer (initializer parameters).
 * - grid: grid to initialize in-place.
 */
typedef void (*GridInitializerFn)(const void* env, Grid* grid);

/* A grid initializer object.
 * - env: opaque parameter storage for initializer logic.
 * - fn:  initializer implementation.
 */
typedef struct
{
    void* env;
    GridInitializerFn fn;
} GridInitializer;

/* A transition rule object.
 * - env: opaque parameter storage for transition logic.
 * - fn:  transition implementation.
 */
typedef struct
{
    void* env;
    TransitionFn fn;
} Transition;


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
typedef struct
{
    Grid currentGenBuf;
    Grid nextGenBuf;
    Transition transition;
    Neighborhood neighborhood;
    int historyIndex;
    Cell history[CELLS * HISTORY_SIZE];
} CellularAutomaton;

CellularAutomaton* ca_alloc_2d(Neighborhood neighborhood);
void ca_free_2d(CellularAutomaton* ca);

/* Grid */
void grid_init_empty(Grid* grid, int dim, int cellCount);
void grid_randomize(Grid* grid, int percentAlive);
void grid_clone(const Grid* srcGrid, Grid* dstGrid);

/* Wolfram 1D */
void wolfram_apply_rule_1d(const Grid* currentGrid, WolframCode ruleCode, Cell* outCells);
void wolfram_run_init(WolframRun* run, WolframCode ruleCode);
void wolfram_run_step(WolframRun* run);

/* Neighborhoods */
void van_neumann_neighborhood(Neighborhood* neighborhood);
void moore_neighborhood(Neighborhood* neighborhood);

/* Core CA helpers / steps */
int ca_count_alive(Cell** neighborCells, size_t neighborCount);
void ca_get_neighbors(CellularAutomaton* automaton, int cellIndex, Cell** out_neighbors);
int ca_apply_bs_rule(const char* b, const char* s, int currentState, int aliveNeighbors);

/* Transitions */
int ca_transition_game_of_life(const void* env, const Cell* currentCell, Cell** neighbors, size_t neighborCount);
int ca_transition_replicator(const void* env, const Cell* currentCell, Cell** neighbors, size_t neighborCount);

/* Allocated generic B/S transition */
Transition* transition_alloc_bs_rule(const char* b, const char* s);
void transition_free_bs_rule(Transition* t);

/* Grid initializer alloc/free */
GridInitializer* grid_initializer_alloc_percent_alive(int percentAlive);
void grid_initializer_free_percent_alive(GridInitializer* gridInitializer);

GridInitializer* grid_initializer_alloc_single_cell_alive(void);
void grid_initializer_free_single_cell_alive(GridInitializer* gridInitializer);

/* Simulation stepping */
void ca_step_forward(CellularAutomaton* ca);
void ca_step_back(CellularAutomaton* automaton);

#ifdef __cplusplus
}
#endif

#endif /* CA_CORE_H */
