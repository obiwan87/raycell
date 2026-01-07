#include "raylib.h"
#include "ca_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

/* Window dimensions in pixels. */
#define WINDOW_WIDTH  800
#define WINDOW_HEIGHT 800

/* Visual size of a single cell in pixels. */
#define CELL_WIDTH   4
#define CELL_HEIGHT  4

/* Background color used for clearing the window. */
#define BACKGROUND ((Color){ 50, 50, 50, 255 })

/* Map a cell state to a render color.
 * NOTE: Currently only handles DEAD/ALIVE. If you add multi-state,
 * update this function (e.g. palette lookup).
 */
static Color raycell_cell_color(int state)
{
    return (state == CELL_STATE_ALIVE)
               ? (Color){180, 80, 10, 255}
               : (Color){50, 50, 50, 255};
}


/* A single automaton configuration slot, combining:
 * - a Transition rule
 * - a GridInitializer for seeding/resetting the grid
 *
 * NOTE: This is used by the "digit-to-select-config" logic in raycell_loop().
 */
typedef struct
{
    Transition* transition;
    GridInitializer* gridInitializer;
} CellularAutomatonConfig;


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
static void raycell_render_1d_ca(WolframRun* run, Camera2D camera, float intervalSeconds, float* timerSeconds)
{
    const float dt = GetFrameTime();

    if (IsKeyDown(KEY_SPACE))
    {
        wolfram_run_step(run);
    }

    if (IsKeyPressed(KEY_R))
    {
        run->current_step = 0;
    }

    *timerSeconds += dt;
    if (*timerSeconds >= intervalSeconds)
    {
        *timerSeconds -= intervalSeconds;
        wolfram_run_step(run);
    }

    const int screen_h = GetScreenHeight();
    const int grid_pixel_h = (run->current_step + 1) * CELL_HEIGHT;

    /* Auto-scroll once the grid becomes taller than a threshold. */
    if (grid_pixel_h * 4 > 3 * screen_h)
    {
        camera.target = (Vector2){0.0f, (float)grid_pixel_h - 3.0f * (float)screen_h / 4.0f};
    }
    else
    {
        camera.target = (Vector2){0.0f, 0.0f};
    }

    BeginDrawing();
    ClearBackground(BACKGROUND);

    BeginMode2D(camera);
    for (int step = 0; step <= run->current_step; step++)
    {
        const Grid* grid = &run->grids[step];
        for (int cell = 0; cell < grid->size; cell++)
        {
            DrawRectangle(cell * CELL_WIDTH,
                          step * CELL_HEIGHT,
                          CELL_WIDTH, CELL_HEIGHT,
                          raycell_cell_color(grid->cells[cell].state));
        }
    }
    EndMode2D();

    EndDrawing();
}


/* Draw the current 2D automaton grid.
 * Assumptions:
 * - The grid is square.
 * - Cells are stored in a 1D array in row-major order.
 * - size == side*side, where side = sqrt(size).
 */
static void raycell_render_2d_ca(CellularAutomaton* automaton)
{
    const int side = (int)sqrt((double)automaton->currentGenBuf.size);

    for (int row = 0; row < side; row++)
    {
        for (int col = 0; col < side; col++)
        {
            const Cell* cell = &automaton->currentGenBuf.cells[row * side + col];
            DrawRectangle(col * CELL_WIDTH,
                          row * CELL_HEIGHT,
                          CELL_WIDTH, CELL_HEIGHT,
                          raycell_cell_color(cell->state));
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
static CellularAutomatonConfig configs[10];
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
static void raycell_loop(void)
{
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
    CellularAutomaton* ca = ca_alloc_2d(neighborhood);

    /* Seed initial grid using config 0 initializer. */
    gofGridInitializer->fn(gofGridInitializer->env, &ca->currentGenBuf);

    
    float interval = 0.1f;
    float timer = 0.0f;
    while (!WindowShouldClose())
    {
        float dt = GetFrameTime(); 
        timer += dt;

        /* Detect digit keys (top row and keypad) to switch configs. */
        int digit = -1;
        for (int i = 0; i <= 9; i++)
        {
            if (IsKeyPressed(KEY_ZERO + i))
            {
                digit = i;
                break;
            }
        }
        for (int i = 0; i <= 9 && digit == -1; i++)
        {
            if (IsKeyPressed(KEY_KP_0 + i))
            {
                digit = i;
                break;
            }
        }

        /* Reset logic: either 'R' or selecting a new config. */
        bool shouldReset = IsKeyPressed(KEY_R);
        if (digit >= 0 && digit < configs_count)
        {
            printf("Digit pressed: %d\n", digit);
            active_config = digit;
            shouldReset = true;
        }

        /* Apply active config to CA. */
        const GridInitializer* activeInitializer = configs[active_config].gridInitializer;
        const Transition* activeTransition = configs[active_config].transition;
        ca->transition = *activeTransition;

        if (shouldReset)
        {
            activeInitializer->fn(activeInitializer->env, &ca->currentGenBuf);
            ca->historyIndex = 0;
        }

        /* Manual stepping. */
        if (IsKeyPressed(KEY_SPACE))
        {
            ca_step_forward(ca);
        }

        if (IsKeyPressed(KEY_LEFT))
        {
            ca_step_back(ca);
        }

        if(timer >= interval) {
            timer -= interval;
            ca_step_forward(ca);
        }

        /* Draw current generation. */
        BeginDrawing();
        ClearBackground(BACKGROUND);
        raycell_render_2d_ca(ca);
        EndDrawing();
    }

    CloseWindow();

    /* Cleanup heap-allocated configuration objects. */
    transition_free_bs_rule(gofTransition);
    transition_free_bs_rule(replicatorTransition);

    grid_initializer_free_percent_alive(gofGridInitializer);
    grid_initializer_free_single_cell_alive(replicatorInitializer);

    ca_free_2d(ca);
}

/* Program entry point. */
int main(void)
{
    /* Seed PRNG for grid randomization and single-cell seed selection. */
    srand(696969);

    raycell_loop();
    return 0;
}
