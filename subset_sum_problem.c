//
// AED, November 2021
//
// Solution of the first practical assignement (subset sum problem)
//
// 103154 - João Fonseca
// 103183 - Diogo Paiva
// 103865 - Jorge Silva
//

#if __STDC_VERSION__ < 199901L
#error "This code must must be compiled in c99 mode or later (-std=c99)" // to handle the unsigned long long data type
#endif
#ifndef STUDENT_H_FILE
#define STUDENT_H_FILE "103183_extra.h"
#endif

/**
 * custom data types
 *
 * the STUDENT_H_FILE defines the following constants and data types
 *
 *   #define min_n       24                   --- the smallest n value we will handle
 *   #define max_n       57                   --- the largest n value we will handle
 *   #define n_sums      20                   --- the number of sums for each n value
 *   #define n_problems  (max_n - min_n + 1)  --- the number of n values
 *
 *   typedef unsigned long long integer_t;    ---  64-bit unsigned integer
 *   typedef struct
 *   {
 *     int n;                                 --- number of elements of the set (for a valid problem, min_n <= n <= max_n)
 *     integer_t p[max_n];                    --- the elements of the set, already sorted in increasing order (only the first n elements are used)
 *     integer_t sums[n_sums];                --- several sums (problem: for each sum find the corresponding subset)
 *   }
 *   subset_sum_problem_data_t;               --- weights p[] and sums for a given value of n
 *
 *   subset_sum_problem_data_t all_subset_sum_problems[n_problems]; --- // the problems
 */

/* ------------------------------------------- Macros ------------------------------------------- */

/*
 * determinar problemas com n pertencente ao intervalo [N_MIN N_MAX]
 */
#define N_MIN 10
#define N_MAX 32

/*
 * 0 - brute force não recursiva
 * 1 - brute force recursiva
 * 2 - brute force clever
 * 3 - horowitz and sahni
 * 4 - schroeppel and shamir
 */
#define FUNC 2

/*
 * utilização correta do programa
 */
#define USAGE "Sintaxe - %s [ficheiro] [w|a]\n"

/*
 * configuração do programa
 */
#define CONFIG "Configuração do programa\n" \
               "    N_MIN ....... %d\n"       \
               "    N_MAX ....... %d\n"       \
               "    FUNC ........ %d\n"       \
               "Problemas do ficheiro %s\n"   \
               "    min_n ....... %d\n"       \
               "    max_n ....... %d\n"       \
               "    n_sums ...... %d\n"       \
               "    n_problems .. %d\n"       \
               "    integer_t ... %d bits\n"

/* ----------------------------------- Inclusão de Ficheiros ------------------------------------ */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "elapsed_time.h"
#include STUDENT_H_FILE

#if FUNC == 0
#include "heap.h"
#endif

/* ------------------------------------ Estruturas de Dados ------------------------------------- */

#if FUNC == 3 || FUNC == 4

typedef struct
{
    integer_t sum;
    unsigned int mask;
} sm_data_t;

#endif

/* ------------------------------------- Funções de Suporte ------------------------------------- */

#if FUNC == 3 || FUNC == 4

/**
 * @brief Gera todas as somas e máscaras dos valores de um subconjunto de p, por um método recursivo.
 *
 * @param p             Conjunto com os valores a somar
 * @param offset        Desvio em relação ao primeiro valor de p
 * @param size          Tamanho do subconjunto
 * @param current_index Índice atual no subconjunto
 * @param partial_sum   Soma atual
 * @param mask          Máscara atual
 * @param r             Array que guarda as somas e as máscaras
 */
void gen_sm_rec(integer_t *p, unsigned int offset, unsigned int size, unsigned int current_index, integer_t partial_sum, unsigned int mask, sm_data_t *r)
{
    /* Ainda não foram percorridos todos os valores do subconjunto */
    if (current_index < size)
    {
        /* p[current_index] não é usado na soma */
        gen_sm_rec(p, offset, size, current_index + 1, partial_sum, mask, r);
        /* p[current_index] é usado na soma */
        gen_sm_rec(p, offset, size, current_index + 1, partial_sum + p[offset + current_index], mask | (1 << current_index), r);
    }
    else
    {
        /* guardar a soma e a máscara no array r */
        r[mask].sum = partial_sum;
        r[mask].mask = mask;
    }
}

int sm_data_cmpfunc(const void *d1, const void *d2)
{
    sm_data_t *data1 = (sm_data_t *)d1;
    sm_data_t *data2 = (sm_data_t *)d2;
    return (data1->sum > data2->sum) - (data1->sum < data2->sum);
}

#endif

/* ----------------------------------- Funções dos Algoritmos ----------------------------------- */

#if FUNC == 0

/**
 * @brief Determina a combinação dos valores de p cujo somatório é desired_sum, por um método brute force não recursivo.
 *
 * @param n             Tamanho de p
 * @param p             Conjunto com os valores a somar
 * @param desired_sum   Soma desejada
 * @param r             Array que guarda a solução
 *
 * @return              1 se foi encontrada uma solução, 0 caso contrário
 */
int brute_force(int n, integer_t p[n], integer_t desired_sum, int r[n])
{
    __uint128_t mask;   // máscara
    integer_t test_sum; // soma de teste

    /* para cada combinação */
    for (mask = 0; mask < (__uint128_t)1 << n; mask++)
    {
        /* determinar a soma dos valores de p */
        test_sum = 0;
        for (int i = 0; i < n; i++)
            if (mask & (1 << i))
                test_sum += p[i];

        /* foi descoberta a soma */
        if (test_sum == desired_sum)
        {
            /* guardar a máscara no array r */
            for (int i = 0; i < n; i++)
                r[i] = (mask & (1 << i)) ? 1 : 0;
            return 1;
        }
    }

    return 0;
}

#elif FUNC == 1

/**
 * @brief Determina a combinação dos valores de p cujo somatório é desired_sum, por um método brute force recursivo.
 *
 * @param n                 Tamanho de p
 * @param p                 Conjunto com os valores a somar
 * @param desired_sum       Soma desejada
 * @param current_index     Índice atual em p
 * @param partial_sum       Soma atual
 * @param mask              Máscara atual
 * @param r                 Array que guarda a solução
 *
 * @return                  1 se foi encontrada uma solução, 0 caso contrário
 */
int brute_force_recursive(int n, integer_t p[n], integer_t desired_sum, int current_index, integer_t partial_sum, unsigned long long mask, int r[n])
{
    /* foi descoberta a soma */
    if (partial_sum == desired_sum)
    {
        /* guardar a máscara no array r */
        for (int i = 0; i < n; i++)
            r[i] = (mask & (1 << i)) ? 1 : 0;
        return 1;
    }

    /* valor de retorno é zero se ainda não tiver sido encontrada uma solução */
    int ret = 0;
    if (current_index < n)
    {
        /* p[current_index] não é usado na soma */
        ret |= brute_force_recursive(n, p, desired_sum, current_index + 1, partial_sum, mask, r);
        /* p[current_index] é usado na soma */
        ret |= brute_force_recursive(n, p, desired_sum, current_index + 1, partial_sum + p[current_index], mask | (1 << current_index), r);
    }
    return ret;
}

#elif FUNC == 2

/**
 * @brief Determina a combinação dos valores de p cujo somatório é desired_sum, por um método brute force recursivo inteligente,
 * que evita recursões extra quando a soma parcial já é superior à soma desejada.
 *
 * @param n                 Tamanho de p
 * @param p                 Conjunto com os valores a somar
 * @param desired_sum       Soma desejada
 * @param current_index     Índice atual em p
 * @param partial_sum       Soma atual
 * @param mask              Máscara atual
 * @param r                 Array que guarda a solução
 *
 * @return                  1 se foi encontrada uma solução, 0 caso contrário
 */
int brute_force_clever(int n, integer_t p[n], integer_t desired_sum, int current_index, integer_t partial_sum, unsigned long long mask, int r[n])
{
    /* foi descoberta a soma */
    if (partial_sum == desired_sum)
    {
        /* guardar a máscara no array r */
        for (int i = 0; i < n; i++)
            r[i] = (mask & (1 << i)) ? 1 : 0;
        return 1;
    }

    /* valor de retorno é zero se ainda não tiver sido encontrada uma solução */
    int ret = 0;
    if (current_index < n && partial_sum < desired_sum)
    {
        /* p[current_index] não é usado na soma */
        ret |= brute_force_clever(n, p, desired_sum, current_index + 1, partial_sum, mask, r);
        /* p[current_index] é usado na soma */
        ret |= brute_force_clever(n, p, desired_sum, current_index + 1, partial_sum + p[current_index], mask | (1 << current_index), r);
    }
    return ret;
}

#elif FUNC == 3

/**
 * @brief Determina a combinação dos valores de p cujo somatório é desired_sum, pelo método de Horowitz e Sahni.
 *
 * @param n                 Tamanho de p
 * @param p                 Conjunto com os valores a somar
 * @param desired_sum       Soma desejada
 * @param r                 Array que guarda a solução
 *
 * @return                  1 se foi encontrada uma solução, 0 caso contrário
 */
int horowitz_and_sahni(int n, integer_t p[n], integer_t desired_sum, int r[n])
{
    /* tamanho das metades de p */
    unsigned int size_p1 = n / 2;
    unsigned int size_p2 = n - size_p1;

    /* tamanho de a e b */
    unsigned long long size_a = 1ull << size_p1; // 2^size_p1
    unsigned long long size_b = 1ull << size_p2; // 2^size_p2

    /* alocar a e b na memória */
    sm_data_t *a = (sm_data_t *)malloc(size_a * sizeof(sm_data_t));
    sm_data_t *b = (sm_data_t *)malloc(size_b * sizeof(sm_data_t));

    /* gerar elementos de a e b */
    gen_sm_rec(p, 0, size_p1, 0, 0, 0, a);
    gen_sm_rec(p, size_p1, size_p2, 0, 0, 0, b);

    /* ordenar (quicksort) os elementos de a e de b */
    qsort(a, size_a, sizeof(sm_data_t), sm_data_cmpfunc);
    qsort(b, size_b, sizeof(sm_data_t), sm_data_cmpfunc);

    unsigned int ret = 0;            // valor de retorno é zero se não for encontrada nenhuma solução
    unsigned long long i = 0;        // indexa os elementos de a
    signed long long j = size_b - 1; // indexa os elementos de b
    integer_t test_sum;              // soma de teste
    unsigned int k;                  // posição do bit da máscara

    /* método de Horowitz e Sahni */
    while (i < size_a && j >= 0)
    {
        test_sum = a[i].sum + b[j].sum;

        if (test_sum < desired_sum)
        {
            i++;
        }
        else if (test_sum > desired_sum)
        {
            j--;
        }
        /* foi descoberta a soma */
        else
        {
            /* guardar a máscara de a[i] em r */
            for (k = 0; k < size_p1; k++)
                r[k] = (a[i].mask & (1 << k)) ? 1 : 0;
            /* guardar a máscara de b[j] em r */
            for (k = 0; k < size_p2; k++)
                r[k + size_p1] = (b[j].mask & (1 << k)) ? 1 : 0;

            ret = 1;
            break;
        }
    }

    /* desalocar a e b da memória */
    free(a);
    free(b);

    return ret;
}

#elif FUNC == 4

/**
 * @brief Determina a combinação dos valores de p cujo somatório é desired_sum, pelo método de Schroeppel e Shamir.
 *
 * @param n                 Tamanho de p
 * @param p                 Conjunto com os valores a somar
 * @param desired_sum       Soma desejada
 * @param r                 Array que guarda a solução
 *
 * @return                  1 se foi encontrada uma solução, 0 caso contrário
 */
int schroeppel_and_shamir(int n, integer_t p[n], integer_t desired_sum, int r[n])
{
    /* tamanho das metades de p */
    unsigned int size_p12 = n / 2;
    unsigned int size_p34 = n - size_p12;

    /* tamanho dos quartos p1, p2, p3 e p4 */
    unsigned int size_p1 = size_p12 / 2;
    unsigned int size_p2 = size_p12 - size_p1;
    unsigned int size_p3 = size_p34 / 2;
    unsigned int size_p4 = size_p34 - size_p3;

    /* tamanho de a1, a2, b1 e b2 */
    unsigned int size_a1 = 1 << size_p1; // 2^size_p1
    unsigned int size_a2 = 1 << size_p2; // 2^size_p2
    unsigned int size_b1 = 1 << size_p3; // 2^size_p3
    unsigned int size_b2 = 1 << size_p4; // 2^size_p4

    /* alocar a1, a2, b1 e b2 na memória */
    sm_data_t *a1 = (sm_data_t *)malloc(size_a1 * sizeof(sm_data_t));
    sm_data_t *a2 = (sm_data_t *)malloc(size_a2 * sizeof(sm_data_t));
    sm_data_t *b1 = (sm_data_t *)malloc(size_b1 * sizeof(sm_data_t));
    sm_data_t *b2 = (sm_data_t *)malloc(size_b2 * sizeof(sm_data_t));

    /* gerar elementos de a1, a2, b1 e b2 */
    gen_sm_rec(p, 0, size_p1, 0, 0, 0, a1);
    gen_sm_rec(p, size_p1, size_p2, 0, 0, 0, a2);
    gen_sm_rec(p, size_p12, size_p3, 0, 0, 0, b1);
    gen_sm_rec(p, size_p12 + size_p3, size_p4, 0, 0, 0, b2);

    /*
     * ordenar (quicksort) os elementos de a2 e b2
     *
     * a1 e b1 serão ordenados na min-heap e na max-heap respetivamente
     */
    qsort(a2, size_a2, sizeof(sm_data_t), sm_data_cmpfunc);
    qsort(b2, size_b2, sizeof(sm_data_t), sm_data_cmpfunc);

    /* criar uma min-heap (a1 e a2) e uma max-heap (b1 e b2) */
    heap_t *min_heap = create(size_a1, 0);
    heap_t *max_heap = create(size_b1, 1);

    /* preencher as heaps com os valores iniciais */
    heap_data_t heap_data;
    for (unsigned int i = 0; i < size_a1; i++)
    {
        heap_data.sum = a1[i].sum + a2[0].sum;
        heap_data.i1 = i;
        heap_data.i2 = 0;
        push(min_heap, heap_data);
    }
    for (unsigned int i = 0; i < size_b1; i++)
    {
        heap_data.sum = b1[i].sum + b2[size_b2 - 1].sum;
        heap_data.i1 = i;
        heap_data.i2 = size_b2 - 1;
        push(max_heap, heap_data);
    }

    unsigned int ret = 0; // valor de retorno é zero se não for encontrada nenhuma solução
    heap_data_t a;        // elemento da min-heap
    heap_data_t b;        // elemento da max-heap
    integer_t test_sum;   // soma de teste
    unsigned int k;       // posição do bit da máscara

    /* método de Schroeppel e Shamir */
    while (min_heap->count > 0 || max_heap->count > 0)
    {
        /* remover elementos no topo das heaps */
        a = pop(min_heap);
        b = pop(max_heap);

        test_sum = a.sum + b.sum;
        if (test_sum < desired_sum)
        {
            /* elemento ainda pode ser inserido na min-heap */
            if (a.i2 + 1 < size_a2)
            {
                a.sum = a1[a.i1].sum + a2[++a.i2].sum;
                push(min_heap, a); // inserir elemento novamente na min-heap
            }
            push(max_heap, b); // inserir elemento novamente na max-heap
        }
        else if (test_sum > desired_sum)
        {
            /* elemento ainda pode ser inserido na max-heap */
            if ((signed int)b.i2 - 1 >= 0)
            {
                b.sum = b1[b.i1].sum + b2[--b.i2].sum;
                push(max_heap, b); // inserir elemento novamente na max-heap
            }
            push(min_heap, a); // inserir elemento novamente na min-heap
        }
        /* foi descoberta a soma */
        else
        {
            /* guardar a máscara de a1[a.i1] em r */
            for (k = 0; k < size_p1; k++)
                r[k] = (a1[a.i1].mask & (1 << k)) ? 1 : 0;
            /* guardar a máscara de a2[a.i2] em r */
            for (k = 0; k < size_p2; k++)
                r[k + size_p1] = (a2[a.i2].mask & (1 << k)) ? 1 : 0;
            /* guardar a máscara de b1[b.i1] em r */
            for (k = 0; k < size_p3; k++)
                r[k + size_p12] = (b1[b.i1].mask & (1 << k)) ? 1 : 0;
            /* guardar a máscara de b2[b.i2] em r */
            for (k = 0; k < size_p4; k++)
                r[k + size_p12 + size_p3] = (b2[b.i2].mask & (1 << k)) ? 1 : 0;

            ret = 1;
            break;
        }
    }

    /* destruir a min-heap e a max-heap */
    destroy(min_heap);
    destroy(max_heap);

    /* desalocar a1, a2, b1 e b2 da memória */
    free(a1);
    free(a2);
    free(b1);
    free(b2);

    return ret;
}

#endif

/* ------------------------------------- Programa Principal ------------------------------------- */

int main(int argc, char *argv[])
{
    /* imprimir utilização incorreta do programa */
    if (argc != 1 && argc != 3)
    {
        fprintf(stderr, USAGE, argv[0]);
        return EXIT_FAILURE;
    }

    int ret;                                                    // valor de retorno das funções
    double t1, t2;                                              // instantes inicial e final
    FILE *out = (argc == 3) ? fopen(argv[1], argv[2]) : stdout; // para onde vai o output do programa

    /* ocorreu um erro ao abrir o ficheiro */
    if (out == NULL)
    {
        fprintf(stderr, "Erro ao abrir o ficheiro de escrita!\n");
        return EXIT_FAILURE;
    }

    /* imprimir configuração do programa */
    fprintf(stderr, CONFIG,
            N_MIN, N_MAX, FUNC,
            STUDENT_H_FILE, min_n, max_n, n_sums, n_problems, 8 * (int)sizeof(integer_t));

    /* imprimir cabeçalho */
    if (argc != 3 || (argc == 3 && strcmp(argv[2], "w") == 0))
        fprintf(out, "%s\t%s\t%s\t%s\n", "n", "retval", "time", "result");

    /*
     * para cada problema
     */
    for (int i = 0; i < n_problems; i++)
    {
        int n = all_subset_sum_problems[i].n;        // número de valores a somar
        integer_t *p = all_subset_sum_problems[i].p; // conjunto com os valores a somar
        int r[n];                                    // resultado

        /* ignorar problemas com n fora do intervalo [N_MIN N_MAX] */
        if (n < N_MIN || n > N_MAX)
            continue;

        /*
         * para cada soma
         */
        for (int j = 0; j < n_sums; j++)
        {
            integer_t desired_sum = all_subset_sum_problems[i].sums[j]; // soma desejada

            /* instante inicial */
            t1 = cpu_time();

            /*
             * executar função
             */
#if FUNC == 0
            ret = brute_force(n, p, desired_sum, r);
#elif FUNC == 1
            ret = brute_force_recursive(n, p, desired_sum, 0, 0, 0, r);
#elif FUNC == 2
            ret = brute_force_clever(n, p, desired_sum, 0, 0, 0, r);
#elif FUNC == 3
            ret = horowitz_and_sahni(n, p, desired_sum, r);
#elif FUNC == 4
            ret = schroeppel_and_shamir(n, p, desired_sum, r);
#endif

            /* instante final */
            t2 = cpu_time();

            /* imprimir n do problema atual, valor de retorno e tempo de execução */
            fprintf(out, "%d\t%d\t%f\t", n, ret, (ret) ? t2 - t1 : NAN);

            /* imprimir resultado */
            if (ret)
                for (int i = 0; i < n; i++)
                    fprintf(out, "%d", r[i]);
            else
                fprintf(out, "Não foi encontrada uma solução!");

            fprintf(out, "\n");
        }
    }

    /* fechar ficheiro */
    if (argc == 3)
        fclose(out);

    return EXIT_SUCCESS;
}
