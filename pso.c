/* Implementação do algoritmo Particle Swarm Optimization (PSO)

   Copyright 2010 Kyriakos Kentzoglanakis

   Este programa é software livre: você pode redistribuí-lo e/ou
   modificá-lo sob os termos da GNU General Public License versão 3
   publicada pela Free Software Foundation.

   Este programa é distribuído na esperança de ser útil, mas SEM
   QUALQUER GARANTIA; sem mesmo a garantia implícita de COMERCIALIZAÇÃO
   ou ADEQUAÇÃO A UM PROPÓSITO ESPECÍFICO. Veja a GNU GPL para mais detalhes.

   Você deve ter recebido uma cópia da GNU GPL junto com este programa.
   Se não, veja <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>   // rand(), malloc(), free()
#include <stdio.h>    // printf()
#include <time.h>     // time()
#include <math.h>     // cos(), pow(), sqrt(), fmod()
#include <float.h>    // DBL_MAX
#include <string.h>   // memmove(), memset()

#include "pso.h"

// ==============================================================
//  Saída "gráfica" no terminal (barra de progresso)
//  - Mantém compatibilidade: não muda a API.
//  - Quando settings->print_every > 0, imprime uma barra atualizando
//    na mesma linha (evita poluir o terminal com muitas linhas).
//  - Funciona em Windows Terminal/PowerShell e também em Linux/macOS.
// ==============================================================
static void pso_print_progress_bar(int step, int steps, double w, double best_err) {
    const int bar_width = 28;
    double frac = 0.0;

    // fração de progresso entre 0 e 1
    if (steps > 0) frac = (double)step / (double)steps;
    if (frac < 0.0) frac = 0.0;
    if (frac > 1.0) frac = 1.0;

    int filled = (int)(frac * bar_width);

    // \r volta para o início da linha; não imprime \n (atualiza a linha atual)
    printf("\r[");
    for (int i = 0; i < bar_width; i++) {
        putchar(i < filled ? '#' : '-');
    }
    printf("] %3d%% | step %d/%d | w=%.2f | best=%.5e",
           (int)(frac * 100.0), step, steps, w, best_err);
    fflush(stdout);
}

// =============================================================
// Macros de números aleatórios
// =============================================================

// gera um double no intervalo (0, 1)
#define RNG_UNIFORM() (rand()/(double)RAND_MAX)

// gera um inteiro no intervalo [0, s)
#define RNG_UNIFORM_INT(s) (rand()%s)

// tipo de função para as diferentes estratégias de vizinhança (informação)
typedef void (*inform_fun_t)(int *comm, double **pos_nb,
                             double **pos_b, double *fit_b,
                             double *gbest, int improved,
                             pso_settings_t *settings);

// tipo de função para as diferentes estratégias de inércia
typedef double (*inertia_fun_t)(int step, pso_settings_t *settings);

//==============================================================
// Calcula tamanho do enxame com base na dimensionalidade
// (heurística comum: 10 + 2*sqrt(dim))
int pso_calc_swarm_size(int dim) {
    int size = 10. + 2. * sqrt(dim);
    return (size > PSO_MAX_SIZE ? PSO_MAX_SIZE : size);
}

//==============================================================
//         ESTRATÉGIAS DE ATUALIZAÇÃO DA INÉRCIA (w)
//==============================================================

// Inércia decrescente linear:
// - começa em w_max e vai caindo até w_min
// - normalmente ajuda a explorar no início e refinar no final
double calc_inertia_lin_dec(int step, pso_settings_t *settings) {

    // fase de decaimento (3/4 do total de steps)
    int dec_stage = 3 * settings->steps / 4;

    if (step <= dec_stage)
        return settings->w_min + (settings->w_max - settings->w_min) *
            (dec_stage - step) / dec_stage;
    else
        return settings->w_min;
}

//==============================================================
//      ESTRATÉGIAS DE VIZINHANÇA (MATRIZ COMM)
//==============================================================

// Vizinhança GLOBAL:
// todas as partículas são informadas pelo melhor global (gbest)
void inform_global(int *comm, double **pos_nb,
                   double **pos_b, double *fit_b,
                   double *gbest, int improved,
                   pso_settings_t *settings)
{
    (void)comm; (void)pos_b; (void)fit_b; (void)improved;
    // todas recebem o mesmo "atrator": gbest
    for (int i=0; i<settings->size; i++)
        memmove((void *)pos_nb[i], (void *)gbest,
                sizeof(double) * settings->dim);
}

// ---------------------------------------------------------------
// Função geral "inform":
// Dada a matriz COMM (conectividade), encontra o melhor "informante"
// para cada partícula e copia a posição pbest do melhor vizinho para pos_nb.
// pos_nb[j] = melhor posição encontrada entre os vizinhos de j
// ---------------------------------------------------------------
void inform(int *comm, double **pos_nb, double **pos_b, double *fit_b,
            int improved, pso_settings_t * settings)
{
    (void)improved;
    int i, j;
    int b_n; // índice do melhor vizinho em termos de fitness

    // para cada partícula j
    for (j=0; j<settings->size; j++) {
        b_n = j; // inicialmente, considera ela mesma como melhor
        // procura qual vizinho (informante) tem menor erro
        for (i=0; i<settings->size; i++)
            // se i informa j e i tem fitness melhor que o melhor atual
            if (comm[i*settings->size + j] && fit_b[i] < fit_b[b_n])
                b_n = i;

        // copia o pbest do melhor vizinho para pos_nb[j]
        memmove((void *)pos_nb[j],
                (void *)pos_b[b_n],
                sizeof(double) * settings->dim);
    }
}

// ============================
// Topologia RING (anel)
// ============================

// Inicializa a matriz COMM para topologia em anel (fixa):
// cada partícula se conecta com ela mesma + vizinho da esquerda + vizinho da direita
void init_comm_ring(int *comm, pso_settings_t * settings) {
    // zera a matriz de conectividade
    memset((void *)comm, 0, sizeof(int)*settings->size*settings->size);

    for (int i=0; i<settings->size; i++) {
        // cada partícula informa a si mesma
        comm[i*settings->size+i] = 1;

        if (i==0) {
            // vizinho à direita
            comm[i*settings->size+i+1] = 1;
            // vizinho à esquerda (último da lista)
            comm[(i+1)*settings->size-1] = 1;
        } else if (i==settings->size-1) {
            // vizinho à direita (primeiro da lista)
            comm[i*settings->size] = 1;
            // vizinho à esquerda
            comm[i*settings->size+i-1] = 1;
        } else {
            // vizinho à direita
            comm[i*settings->size+i+1] = 1;
            // vizinho à esquerda
            comm[i*settings->size+i-1] = 1;
        }
    }
}

void inform_ring(int *comm, double **pos_nb,
                 double **pos_b, double *fit_b,
                 double *gbest, int improved,
                 pso_settings_t * settings)
{
    (void)gbest;
    // atualiza pos_nb usando a matriz COMM do anel
    inform(comm, pos_nb, pos_b, fit_b, improved, settings);
}

// ============================
// Topologia RANDOM (aleatória)
// ============================

// Inicializa COMM de forma aleatória:
// em média, cada partícula escolhe nhood_size informantes
void init_comm_random(int *comm, pso_settings_t * settings) {
    // zera a matriz
    memset((void *)comm, 0, sizeof(int)*settings->size*settings->size);

    for (int i=0; i<settings->size; i++) {
        // cada partícula informa a si mesma
        comm[i*settings->size + i] = 1;

        // escolhe informantes aleatórios
        for (int k=0; k<settings->nhood_size; k++) {
            int j = RNG_UNIFORM_INT(settings->size);
            // partícula i informa partícula j
            comm[i*settings->size + j] = 1;
        }
    }
}

void inform_random(int *comm, double **pos_nb,
                   double **pos_b, double *fit_b,
                   double *gbest, int improved,
                   pso_settings_t * settings)
{
    (void)gbest;

    // Se não houve melhora, muda (regenera) a vizinhança aleatória
    if (!improved)
        init_comm_random(comm, settings);

    inform(comm, pos_nb, pos_b, fit_b, improved, settings);
}

//==============================================================
// Cria configurações padrão do PSO
pso_settings_t *pso_settings_new(int dim, double range_lo, double range_hi) {
    pso_settings_t *settings = (pso_settings_t *)malloc(sizeof(pso_settings_t));
    if (settings == NULL) { return NULL; }

    // valores padrão
    settings->dim = dim;
    settings->goal = 1e-5;

    // aloca e inicializa limites por dimensão
    settings->range_lo = (double *)malloc(settings->dim * sizeof(double));
    if (settings->range_lo == NULL) { free(settings); return NULL; }

    settings->range_hi = (double *)malloc(settings->dim * sizeof(double));
    if (settings->range_hi == NULL) { free(settings); free(settings->range_lo); return NULL; }

    for (int i=0; i<settings->dim; i++) {
        settings->range_lo[i] = range_lo;
        settings->range_hi[i] = range_hi;
    }

    // defaults clássicos
    settings->size = pso_calc_swarm_size(settings->dim);
    settings->print_every = 1000;
    settings->steps = 100000;
    settings->c1 = 1.496;
    settings->c2 = 1.496;
    settings->w_max = PSO_INERTIA;
    settings->w_min = 0.3;

    settings->clamp_pos = 1;
    settings->nhood_strategy = PSO_NHOOD_RING;
    settings->nhood_size = 5;
    settings->w_strategy = PSO_W_LIN_DEC;

    return settings;
}

// Libera configurações do PSO
void pso_settings_free(pso_settings_t *settings) {
    free(settings->range_lo);
    free(settings->range_hi);
    free(settings);
}

//==============================================================
// Funções auxiliares: criação/liberação de matrizes
double **pso_matrix_new(int size, int dim) {
    double **m = (double **)malloc(size * sizeof(double *));
    for (int i=0; i<size; i++) {
        m[i] = (double *)malloc(dim * sizeof(double));
    }
    return m;
}

void pso_matrix_free(double **m, int size) {
    for (int i=0; i<size; i++) {
        free(m[i]);
    }
    free(m);
}

//==============================================================
//                 ALGORITMO PSO (PRINCIPAL)
//==============================================================
void pso_solve(pso_obj_fun_t obj_fun, void *obj_fun_params,
               pso_result_t *solution, pso_settings_t *settings)
{
    // ---------------------------
    // Estruturas das partículas
    // ---------------------------
    // pos   : posições atuais
    // vel   : velocidades atuais
    // pos_b : melhor posição (pbest) de cada partícula
    double **pos   = pso_matrix_new(settings->size, settings->dim);
    double **vel   = pso_matrix_new(settings->size, settings->dim);
    double **pos_b = pso_matrix_new(settings->size, settings->dim);

    // fit   : fitness (erro) atual de cada partícula
    // fit_b : melhor fitness (erro) de cada partícula (pbest)
    double *fit   = (double *)malloc(settings->size * sizeof(double));
    double *fit_b = (double *)malloc(settings->size * sizeof(double));

    // pos_nb : melhor posição informada (melhor dos vizinhos) para cada partícula
    double **pos_nb = pso_matrix_new(settings->size, settings->dim);

    // comm : matriz de conectividade (quem informa quem)
    int *comm = (int *)malloc(settings->size * settings->size * sizeof(int));

    // improved indica se o gbest melhorou na última iteração
    int improved = 0;

    int i, d, step;
    double a, b;       // usados na inicialização (posição/velocidade)
    double rho1, rho2; // coeficientes aleatórios
    double w = PSO_INERTIA; // inércia atual

    inform_fun_t  inform_fun = NULL;     // função de vizinhança
    inertia_fun_t calc_inertia_fun = NULL; // função de inércia

    // semente aleatória
    srand(time(NULL));

    // ---------------------------
    // Escolhe a estratégia de vizinhança
    // ---------------------------
    switch (settings->nhood_strategy) {
        case PSO_NHOOD_GLOBAL:
            inform_fun = inform_global;
            break;
        case PSO_NHOOD_RING:
            init_comm_ring(comm, settings);
            inform_fun = inform_ring;
            break;
        case PSO_NHOOD_RANDOM:
            init_comm_random(comm, settings);
            inform_fun = inform_random;
            break;
        default:
            inform_fun = inform_global;
            break;
    }

    // ---------------------------
    // Escolhe a estratégia de inércia
    // ---------------------------
    switch (settings->w_strategy) {
        case PSO_W_LIN_DEC:
            calc_inertia_fun = calc_inertia_lin_dec;
            break;
        default:
            // se não definido, fica como constante (w = PSO_INERTIA)
            calc_inertia_fun = NULL;
            break;
    }

    // Inicializa solução (gbest)
    solution->error = DBL_MAX;

    // ---------------------------
    // Inicialização do enxame
    // ---------------------------
    for (i=0; i<settings->size; i++) {
        for (d=0; d<settings->dim; d++) {
            // sorteia dois valores no intervalo [range_lo, range_hi]
            a = settings->range_lo[d] + (settings->range_hi[d] - settings->range_lo[d]) * RNG_UNIFORM();
            b = settings->range_lo[d] + (settings->range_hi[d] - settings->range_lo[d]) * RNG_UNIFORM();

            // posição inicial
            pos[i][d] = a;

            // pbest começa igual à posição inicial
            pos_b[i][d] = a;

            // velocidade inicial (diferença entre dois pontos / 2)
            vel[i][d] = (a-b) / 2.0;
        }

        // calcula fitness inicial
        fit[i] = obj_fun(pos[i], settings->dim, obj_fun_params);
        fit_b[i] = fit[i];

        // atualiza gbest se necessário
        if (fit[i] < solution->error) {
            solution->error = fit[i];
            memmove((void *)solution->gbest, (void *)pos[i],
                    sizeof(double) * settings->dim);
        }
    }

    // ---------------------------
    // Loop principal do PSO
    // ---------------------------
    int progress_used = 0; // indica se usamos a barra (para finalizar com \n)

    for (step=0; step<settings->steps; step++) {
        // registra o passo atual (caso seja usado fora)
        settings->step = step;

        // atualiza inércia (se houver estratégia definida)
        if (calc_inertia_fun != NULL) {
            w = calc_inertia_fun(step, settings);
        }

        // critério de parada: atingiu o objetivo (goal)
        if (solution->error <= settings->goal) {
            if (settings->print_every) {
                if (progress_used) printf("\n");
                printf("Goal achieved @ step %d (error=%.3e) :-)\n", step, solution->error);
            }
            break;
        }

        // encontra o melhor vizinho (pos_nb) para cada partícula
        inform_fun(comm, pos_nb, pos_b, fit_b, solution->gbest, improved, settings);
        improved = 0; // reseta flag

        // atualiza todas as partículas
        for (i=0; i<settings->size; i++) {
            for (d=0; d<settings->dim; d++) {
                // coeficientes estocásticos
                rho1 = settings->c1 * RNG_UNIFORM();
                rho2 = settings->c2 * RNG_UNIFORM();

                // atualização de velocidade (fórmula do PSO)
                vel[i][d] = w * vel[i][d]
                    + rho1 * (pos_b[i][d] - pos[i][d])
                    + rho2 * (pos_nb[i][d] - pos[i][d]);

                // atualização de posição
                pos[i][d] += vel[i][d];

                // tratamento de limites
                if (settings->clamp_pos) {
                    // CLAMP: trava nas bordas e zera velocidade na dimensão
                    if (pos[i][d] < settings->range_lo[d]) {
                        pos[i][d] = settings->range_lo[d];
                        vel[i][d] = 0;
                    } else if (pos[i][d] > settings->range_hi[d]) {
                        pos[i][d] = settings->range_hi[d];
                        vel[i][d] = 0;
                    }
                } else {
                    // PERIÓDICO: “dá a volta” quando ultrapassa limites
                    if (pos[i][d] < settings->range_lo[d]) {
                        pos[i][d] = settings->range_hi[d] - fmod(settings->range_lo[d] - pos[i][d],
                                                                 settings->range_hi[d] - settings->range_lo[d]);
                        vel[i][d] = 0;
                    } else if (pos[i][d] > settings->range_hi[d]) {
                        pos[i][d] = settings->range_lo[d] + fmod(pos[i][d] - settings->range_hi[d],
                                                                 settings->range_hi[d] - settings->range_lo[d]);
                        vel[i][d] = 0;
                    }
                }
            }

            // avalia fitness na nova posição
            fit[i] = obj_fun(pos[i], settings->dim, obj_fun_params);

            // atualiza pbest (melhor pessoal)
            if (fit[i] < fit_b[i]) {
                fit_b[i] = fit[i];
                memmove((void *)pos_b[i], (void *)pos[i],
                        sizeof(double) * settings->dim);
            }

            // atualiza gbest (melhor global)
            if (fit[i] < solution->error) {
                improved = 1;
                solution->error = fit[i];
                memmove((void *)solution->gbest, (void *)pos[i],
                        sizeof(double) * settings->dim);
            }
        }

        // imprime progresso a cada N passos
        if (settings->print_every && (step % settings->print_every == 0)) {
            pso_print_progress_bar(step, settings->steps, w, solution->error);
            progress_used = 1;
        }
    }

    // garante que o prompt não fique "colado" na barra
    if (progress_used) printf("\n");

    // ---------------------------
    // Libera memória
    // ---------------------------
    pso_matrix_free(pos, settings->size);
    pso_matrix_free(vel, settings->size);
    pso_matrix_free(pos_b, settings->size);
    pso_matrix_free(pos_nb, settings->size);
    free(comm);
    free(fit);
    free(fit_b);
}
