/* Implementação do algoritmo Particle Swarm Optimization (PSO)

   Copyright 2010 Kyriakos Kentzoglanakis

   Este programa é software livre: você pode redistribuí-lo e/ou
   modificá-lo sob os termos da GNU General Public License versão 3,
   publicada pela Free Software Foundation.

   Este programa é distribuído na esperança de ser útil, mas SEM
   QUALQUER GARANTIA; sem mesmo a garantia implícita de COMERCIALIZAÇÃO
   ou ADEQUAÇÃO A UM PROPÓSITO ESPECÍFICO. Veja a GNU GPL para mais detalhes.

   Você deve ter recebido uma cópia da GNU GPL junto com este programa.
   Se não, veja <http://www.gnu.org/licenses/>.
*/

#ifndef PSO_H_
#define PSO_H_

// =============================================================
//                     CONSTANTES GERAIS
// =============================================================

// Tamanho máximo permitido para o enxame (swarm)
// (limite de segurança para evitar alocações muito grandes)
#define PSO_MAX_SIZE 100

// Valor padrão da inércia (w) sugerido na literatura
// (referência: Clerc 2002 / constriction factor)
#define PSO_INERTIA 0.7298


// =============================================================
//                 ESQUEMAS DE VIZINHANÇA (NHOOD)
// =============================================================

// 0) Topologia GLOBAL:
// todas as partículas usam o melhor global (gbest) como referência.
// Converge rápido, mas pode ter convergência prematura.
#define PSO_NHOOD_GLOBAL 0

// 1) Topologia RING (anel):
// cada partícula se comunica com poucos vizinhos (esquerda/direita).
// Ajuda a manter diversidade e evita travar cedo.
#define PSO_NHOOD_RING 1

// 2) Topologia RANDOM (aleatória):
// os informantes mudam ao longo do tempo (dependendo de melhora).
// Pode aumentar a exploração do espaço.
// Referência: http://clerc.maurice.free.fr/pso/random_topology.pdf
#define PSO_NHOOD_RANDOM 2


// =============================================================
//           ESTRATÉGIAS DE ATUALIZAÇÃO DA INÉRCIA (w)
// =============================================================

// 0) Inércia constante (w fixo)
#define PSO_W_CONST 0

// 1) Inércia decrescente linear (w cai de w_max para w_min)
#define PSO_W_LIN_DEC 1


// =============================================================
//              ESTRUTURA DE RESULTADO DO PSO
// =============================================================
// OBS: Esta estrutura deve ser preparada pelo usuário antes de chamar pso_solve().
// O ponteiro gbest deve ser alocado com DIM elementos.
typedef struct {

    // Melhor erro (valor mínimo da função objetivo encontrado)
    double error;

    // Melhor posição global encontrada (gbest)
    // Deve ter exatamente DIM elementos
    double *gbest;

} pso_result_t;


// =============================================================
//              TIPO DA FUNÇÃO OBJETIVO (OBJ FUN)
// =============================================================
// A função objetivo recebe:
// - ponteiro para vetor de posição (double *x)
// - dimensão do problema (int dim)
// - ponteiro genérico para parâmetros extras (void *params)
// E retorna um double: o valor de erro/fitness (quanto menor, melhor).
typedef double (*pso_obj_fun_t)(double *, int, void *);


// =============================================================
//                ESTRUTURA DE CONFIGURAÇÃO DO PSO
// =============================================================
typedef struct {

    // Dimensão do problema (número de variáveis)
    int dim;

    // Limites inferiores e superiores por dimensão (arrays de tamanho DIM)
    double *range_lo; // limite inferior
    double *range_hi; // limite superior

    // Objetivo de parada (se erro <= goal, o PSO para)
    double goal;

    // Tamanho do enxame (número de partículas)
    int size;

    // Frequência de impressão (a cada N passos)
    // Se 0, não imprime nada durante a execução
    int print_every;

    // Número máximo de iterações
    int steps;

    // Passo atual (a biblioteca atualiza isso internamente)
    int step;

    // Coeficientes do PSO:
    // c1 = componente cognitiva (tendência a voltar ao pbest)
    // c2 = componente social (tendência a ir ao gbest ou melhor vizinho)
    double c1;
    double c2;

    // Parâmetros de inércia:
    // w_max e w_min são usados quando a estratégia é linear decrescente
    double w_max;
    double w_min;

    // Controle de limites:
    // clamp_pos = 1 → trava nas bordas (CLAMP) e zera velocidade na dimensão
    // clamp_pos = 0 → condição periódica (quando sai, "dá a volta")
    int clamp_pos;

    // Estratégia de vizinhança:
    // PSO_NHOOD_GLOBAL, PSO_NHOOD_RING ou PSO_NHOOD_RANDOM
    int nhood_strategy;

    // Tamanho médio da vizinhança (usado em RANDOM / e pode ser usado em outras)
    int nhood_size;

    // Estratégia de inércia:
    // PSO_W_CONST ou PSO_W_LIN_DEC
    int w_strategy;

} pso_settings_t;


// =============================================================
//                  FUNÇÕES PÚBLICAS DA BIBLIOTECA
// =============================================================

// Cria e inicializa a estrutura de configurações do PSO
// dim: dimensão do problema
// range_lo / range_hi: limites inferiores/superiores (mesmo valor para todas as dimensões)
pso_settings_t *pso_settings_new(int dim, double range_lo, double range_hi);

// Libera a memória da estrutura de configurações
void pso_settings_free(pso_settings_t *settings);

// Retorna um tamanho de enxame sugerido com base na dimensão
int pso_calc_swarm_size(int dim);

// Executa o PSO para minimizar a função objetivo obj_fun
// - obj_fun: função a ser minimizada
// - obj_fun_params: parâmetros extras (pode ser NULL)
// - solution: onde o PSO escreve o melhor resultado (gbest e error)
// - settings: parâmetros do PSO
void pso_solve(pso_obj_fun_t obj_fun, void *obj_fun_params,
               pso_result_t *solution, pso_settings_t *settings);

#endif // PSO_H_
