#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "pso.h"

#ifdef _WIN32
#include <windows.h>
#endif

// ============================
//   CORES ANSI (Terminal)
// ============================
#define C_RESET     "\x1b[0m"
#define C_DIM       "\x1b[2m"
#define C_BOLD      "\x1b[1m"
#define C_TITLE     "\x1b[1;34m"  // azul escuro (negrito)
#define C_BLUE_DARK "\x1b[34m"    // azul escuro
#define C_GREEN     "\x1b[32m"
#define C_YELLOW    "\x1b[33m"
#define C_RED       "\x1b[31m"
#define C_WHITE     "\x1b[37m"

static void enable_ansi(void) {
#ifdef _WIN32
    HANDLE hOut = GetStdHandle(STD_OUTPUT_HANDLE);
    if (hOut == INVALID_HANDLE_VALUE) return;
    DWORD mode = 0;
    if (!GetConsoleMode(hOut, &mode)) return;
    mode |= 0x0004; // ENABLE_VIRTUAL_TERMINAL_PROCESSING
    SetConsoleMode(hOut, mode);
#endif
}

static void clear_screen(void) {
#ifdef _WIN32
    system("cls");
#else
    system("clear");
#endif
}

static void pause_enter(void) {
    int c;
    printf(C_DIM "\n[ENTER] para continuar..." C_RESET);
    while ((c=getchar())!='\n' && c!=EOF) {}
    getchar();
}

// ============================
//  UTIL: CAIXAS ALINHADAS
// ============================
static const int W = 62; // largura total da caixa (inclui + e |)

static void box_border(char left, char mid, char right, char fill) {
    putchar(left);
    for (int i = 0; i < W - 2; i++) putchar(fill);
    putchar(right);
    putchar('\n');
}

static void box_line_plain(const char *text) {
    // imprime: | <text> ...padding... |
    int inner = W - 4; // espaço útil entre "| " e " |"
    int len = (int)strlen(text);
    if (len > inner) len = inner;

    printf("| ");
    fwrite(text, 1, len, stdout);
    for (int i = 0; i < inner - len; i++) putchar(' ');
    printf(" |\n");
}

static void box_line_color(const char *color, const char *text) {
    int inner = W - 4;
    int len = (int)strlen(text);
    if (len > inner) len = inner;

    printf("| ");
    printf("%s", color);
    fwrite(text, 1, len, stdout);
    printf("%s", C_RESET);

    for (int i = 0; i < inner - len; i++) putchar(' ');
    printf(" |\n");
}

static void box_blank(void) {
    box_line_plain("");
}

// ============================
//      FUNÇÕES OBJETIVO
// ============================
double pso_sphere(double *x, int dim, void *p) {
    (void)p;
    double s=0.0;
    for(int i=0;i<dim;i++) s += x[i]*x[i];
    return s;
}

double pso_rosenbrock(double *x, int dim, void *p) {
    (void)p;
    if (dim < 2) return 1e9; // evita caso degenerado
    double s=0.0;
    for(int i=0;i<dim-1;i++){
        double a = x[i+1] - x[i]*x[i];
        double b = 1.0 - x[i];
        s += 100.0*a*a + b*b;
    }
    return s;
}

double pso_griewank(double *x, int dim, void *p) {
    (void)p;
    double sum=0.0, prod=1.0;
    for(int i=0;i<dim;i++){
        sum += x[i]*x[i];
        prod *= cos(x[i]/sqrt(i+1.0));
    }
    return sum/4000.0 - prod + 1.0;
}

double pso_rastrigin(double *x, int dim, void *p) {
    (void)p;
    double s = 10.0*dim;
    for(int i=0;i<dim;i++)
        s += x[i]*x[i] - 10.0*cos(2.0*M_PI*x[i]);
    return s;
}

double pso_ackley(double *x, int dim, void *p) {
    (void)p;
    double a=20.0, b=0.2, c=2.0*M_PI;
    double s1=0.0, s2=0.0;
    for(int i=0;i<dim;i++){
        s1 += x[i]*x[i];
        s2 += cos(c*x[i]);
    }
    return -a*exp(-b*sqrt(s1/dim)) - exp(s2/dim) + a + exp(1.0);
}

// ============================
//   UI: HEADER / MENU / CARD
// ============================
static void header_box(void) {
    box_border('+','-','+','-');
    box_blank();
    box_line_color(C_TITLE, "Particle Swarm Optimization (PSO)");
    box_line_color(C_BLUE_DARK, "Demonstracao em C");
    box_blank();
    box_border('+','-','+','-');
    putchar('\n');
}

static void menu_box(void) {
    box_border('+','-','+','-');
    box_line_plain("MENU PRINCIPAL");
    box_border('+','-','+','-');

    // opções (sem quebrar borda: texto simples dentro)
    box_line_plain("1 - Sphere     (Esfera - simples)");
    box_line_plain("2 - Rosenbrock (Vale estreito)");
    box_line_plain("3 - Griewank   (Muitos minimos locais)");
    box_line_plain("4 - Rastrigin  (Altamente multimodal)");
    box_line_plain("5 - Ackley     (Multimodal complexa)");
    box_blank();
    box_line_plain("9 - Configurar parametros");
    box_line_plain("0 - Sair");

    box_border('+','-','+','-');
}

static void config_card(int dim, int particles, int steps, double goal,
                        int topo, int wstrat, int clamp, double c1, double c2, int print_every) {
    char buf[128];

    putchar('\n');
    box_border('+','-','+','-');
    box_line_plain("CONFIGURACAO ATUAL");
    box_border('+','-','+','-');

    snprintf(buf, sizeof(buf), "Dimensao   : %d", dim);                 box_line_plain(buf);
    snprintf(buf, sizeof(buf), "Particulas : %d", particles);          box_line_plain(buf);
    snprintf(buf, sizeof(buf), "Steps      : %d", steps);              box_line_plain(buf);
    snprintf(buf, sizeof(buf), "Goal       : %.1e", goal);             box_line_plain(buf);

    snprintf(buf, sizeof(buf), "Topologia  : %s", topo==0?"GLOBAL":"RING"); box_line_plain(buf);
    snprintf(buf, sizeof(buf), "Inercia    : %s", wstrat==0?"CONST":"LIN_DEC"); box_line_plain(buf);
    snprintf(buf, sizeof(buf), "Limites    : %s", clamp? "CLAMP":"PERIODICO"); box_line_plain(buf);

    snprintf(buf, sizeof(buf), "c1/c2      : %.3f / %.3f", c1, c2);    box_line_plain(buf);
    snprintf(buf, sizeof(buf), "Print      : a cada %d passos", print_every); box_line_plain(buf);

    box_border('+','-','+','-');
}

// ============================
// Leitura segura (ENTER = padrão)
// ============================
static int read_int(const char *prompt, int minv, int maxv, int defv) {
    char buf[64];
    int v;
    for (;;) {
        printf("%s%s (min=%d max=%d, padrao=%d)%s: ",
               C_WHITE, prompt, minv, maxv, defv, C_RESET);
        if (!fgets(buf, sizeof(buf), stdin)) return defv;
        if (buf[0] == '\n') return defv;
        if (sscanf(buf, "%d", &v) == 1 && v >= minv && v <= maxv) return v;
        printf(C_RED "Entrada invalida. Tente novamente.\n" C_RESET);
    }
}

static double read_double(const char *prompt, double minv, double maxv, double defv) {
    char buf[64];
    double v;
    for (;;) {
        printf("%s%s (min=%.2g max=%.2g, padrao=%.2g)%s: ",
               C_WHITE, prompt, minv, maxv, defv, C_RESET);
        if (!fgets(buf, sizeof(buf), stdin)) return defv;
        if (buf[0] == '\n') return defv;
        if (sscanf(buf, "%lf", &v) == 1 && v >= minv && v <= maxv) return v;
        printf(C_RED "Entrada invalida. Tente novamente.\n" C_RESET);
    }
}

// ============================
//            MAIN
// ============================
int main(void) {
    enable_ansi();

    // Configuração padrão (boa para apresentação)
    int dim = 30;
    int particles = 30;
    int steps = 10000;
    double goal = 1e-5;
    int topo = 1;          // 0=GLOBAL  1=RING
    int wstrat = 1;        // 0=CONST   1=LIN_DEC
    int clamp = 1;         // 1=CLAMP   0=PERIODICO
    double c1 = 1.496;
    double c2 = 1.496;
    int print_every = 1000;

    for (;;) {
        clear_screen();
        header_box();
        menu_box();
        config_card(dim, particles, steps, goal, topo, wstrat, clamp, c1, c2, print_every);

        printf("\nOpcao: ");
        char in[32];
        if (!fgets(in, sizeof(in), stdin)) break;
        int op = -1;
        sscanf(in, "%d", &op);

        if (op == 0) break;

        if (op == 9) {
            clear_screen();
            header_box();
            printf(C_BOLD "CONFIGURAR PARAMETROS (ENTER = padrao)\n\n" C_RESET);

            dim       = read_int("Dimensao (dim)", 2, 200, dim);
            particles = read_int("Particulas (swarm size)", 10, 200, particles);
            steps     = read_int("Maximo de iteracoes (steps)", 100, 200000, steps);
            goal      = read_double("Goal (parada)", 1e-12, 1e3, goal);

            printf("\nTopologia (0=GLOBAL, 1=RING)\n");
            topo = read_int("Escolha", 0, 1, topo);

            printf("\nInercia (0=CONST, 1=LIN_DEC)\n");
            wstrat = read_int("Escolha", 0, 1, wstrat);

            printf("\nLimites (1=CLAMP, 0=PERIODICO)\n");
            clamp = read_int("Escolha", 0, 1, clamp);

            c1 = read_double("Coeficiente cognitivo (c1)", 0.1, 4.0, c1);
            c2 = read_double("Coeficiente social (c2)",    0.1, 4.0, c2);

            print_every = read_int("Imprimir a cada N passos", 0, 50000, print_every);

            printf(C_GREEN "\nConfiguracao atualizada!\n" C_RESET);
            pause_enter();
            continue;
        }

        // Seleção da função
        pso_obj_fun_t obj = NULL;
        const char *fname = NULL;
        double lo = -100, hi = 100;

        switch (op) {
            case 1: obj = pso_sphere;    fname = "Sphere (Esfera)";     lo=-100;   hi=100;   break;
            case 2: obj = pso_rosenbrock;fname = "Rosenbrock";          lo=-2.048; hi=2.048; break;
            case 3: obj = pso_griewank;  fname = "Griewank";            lo=-600;   hi=600;   break;
            case 4: obj = pso_rastrigin; fname = "Rastrigin";           lo=-5.12;  hi=5.12;  break;
            case 5: obj = pso_ackley;    fname = "Ackley";              lo=-32.0;  hi=32.0;  break;
            default:
                printf(C_RED "\nOpcao invalida.\n" C_RESET);
                pause_enter();
                continue;
        }

        // Preparar settings
        pso_settings_t *settings = pso_settings_new(dim, lo, hi);
        settings->size = particles;
        settings->steps = steps;
        settings->goal = goal;
        settings->c1 = c1;
        settings->c2 = c2;
        settings->print_every = print_every;

        settings->nhood_strategy = (topo==0) ? PSO_NHOOD_GLOBAL : PSO_NHOOD_RING;
        settings->nhood_size = (topo==0) ? particles : (particles < 10 ? particles : 10);

        settings->w_strategy = (wstrat==0) ? PSO_W_CONST : PSO_W_LIN_DEC;
        if (settings->w_strategy == PSO_W_LIN_DEC) {
            settings->w_max = 0.9;
            settings->w_min = 0.4;
        }

        settings->clamp_pos = clamp ? 1 : 0;

        // Resultado
        pso_result_t result;
        result.gbest = (double*)malloc(dim*sizeof(double));
        if (!result.gbest) {
            printf(C_RED "Falha ao alocar memoria.\n" C_RESET);
            pso_settings_free(settings);
            pause_enter();
            continue;
        }

        clear_screen();
        header_box();

        printf(C_BOLD "Funcao: " C_GREEN "%s" C_RESET "\n", fname);
        printf("dim=%d | particulas=%d | steps=%d | goal=%.1e\n", dim, particles, steps, goal);
        printf("topologia=%s | inercia=%s | limites=%s | c1/c2=%.3f/%.3f\n\n",
               topo==0?"GLOBAL":"RING",
               wstrat==0?"CONST":"LIN_DEC",
               clamp? "CLAMP":"PERIODICO",
               c1, c2);

        // Rodar PSO
        pso_solve(obj, NULL, &result, settings);

        // Resultado final (em caixa alinhada)
        box_border('+','-','+','-');
        box_line_plain("RESULTADO");
        box_border('+','-','+','-');

        char linebuf[256];
        snprintf(linebuf, sizeof(linebuf), "Best error : %.12e", result.error);
        box_line_plain(linebuf);

        int show = dim < 10 ? dim : 10;
        // monta gbest curto
        char gbuf[512];
        int off = snprintf(gbuf, sizeof(gbuf), "gbest[%d] : [", show);
        for (int i=0;i<show && off < (int)sizeof(gbuf)-32;i++) {
            off += snprintf(gbuf+off, sizeof(gbuf)-off, "%s%.6f", i?", ":"", result.gbest[i]);
        }
        if (dim > show) snprintf(gbuf+off, sizeof(gbuf)-off, ", ...]");
        else snprintf(gbuf+off, sizeof(gbuf)-off, "]");

        box_line_plain(gbuf);
        box_border('+','-','+','-');

        free(result.gbest);
        pso_settings_free(settings);

        pause_enter();
    }

    return 0;
}
