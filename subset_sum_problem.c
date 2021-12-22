//
// AED, November 2021
//
// Solution of the first practical assignement (subset sum problem)
//
// Place your student numbers and names here
//

#if __STDC_VERSION__ < 199901L
#error "This code must must be compiled in c99 mode or later (-std=c99)" // to handle the unsigned long long data type
#endif
#ifndef STUDENT_H_FILE
#define STUDENT_H_FILE "000000.h"
#endif

//
// include files
//
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "elapsed_time.h"
#include STUDENT_H_FILE

//
// custom data types
//
// the STUDENT_H_FILE defines the following constants and data types
//
//   #define min_n       24                   --- the smallest n value we will handle
//   #define max_n       57                   --- the largest n value we will handle
//   #define n_sums      20                   --- the number of sums for each n value
//   #define n_problems  (max_n - min_n + 1)  --- the number of n values
//
//   typedef unsigned long long integer_t;    ---  64-bit unsigned integer
//   typedef struct
//   {
//     int n;                                 --- number of elements of the set (for a valid problem, min_n <= n <= max_n)
//     integer_t p[max_n];                    --- the elements of the set, already sorted in increasing order (only the first n elements are used)
//     integer_t sums[n_sums];                --- several sums (problem: for each sum find the corresponding subset)
//   }
//   subset_sum_problem_data_t;               --- weights p[] and sums for a given value of n
//
//   subset_sum_problem_data_t all_subset_sum_problems[n_problems]; --- // the problems
//

//
// place your code here
//
// possible function prototype for a recursive brute-force function:
// int brute_force(int n,integer_t p[n],int desired_sum,int current_index,integer_t partial_sum)
// it sould return 1 when the solution is found and 0 otherwise
// note, however, that you may get a faster function by reducing the number of function arguments (maybe a single pointer to a struct?)
//

/* ------------------------------------------- Macros ------------------------------------------- */

/*
 * 1 - imprime na consola o valor de retorno, o tempo de execução e o resultado de cada soma
 * 0 - imprime na consola os tempos de execução das somas de um problema n, na mesma linha e separados por '\t'
 */
#ifndef DEBUG
#define DEBUG 1
#endif

/*
 * Determinar até ao problema com N_LIMIT valores a somar
 */
#ifndef N_LIMIT
#define N_LIMIT 10
#endif

/*
 * 0 - brute force não recursiva
 * 1 - brute force recursiva
 * 2 - brute force clever
 * 3 - horowitz and sahni
 * 4 - schroeppel and shamir
 */
#ifndef FUNC
#define FUNC 3
#endif

/* ------------------------------------ Estruturas de Dados ------------------------------------- */

typedef struct
{
    integer_t sum;
    unsigned long mask;
} hs_data_t;

/* ------------------------------------- Funções de Suporte ------------------------------------- */

int hs_data_cmpfunc(const void *d1, const void *d2)
{
    return ( (hs_data_t*)d1.sum - (hs_data_t*)d2.sum);
}

/* ----------------------------------- Funções dos Algoritmos ----------------------------------- */

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
    unsigned long mask; // máscara
    integer_t test_sum; // soma de teste

    /* para cada combinação */
    for (mask = 0; mask < (1 << n); mask++)
    {
        /* determinar a soma dos valores de p */
        test_sum = 0;
        for (int i = 0; i < n; i++)
            if ((mask & (1 << (n-i-1))) != 0)
                test_sum += p[i];

        /* foi descoberta a soma */
        if (test_sum == desired_sum)
        {
            /* guardar a máscara no array b */
            for (int i = 0; i < n; i++)
                r[i] = ((mask & (1 << (n-i-1))) == 0) ? 0 : 1;
            return 1;
        }
    }

    return 0;
}

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
 *
 */
int brute_force_recursive(int n, integer_t p[n], integer_t desired_sum, int current_index, integer_t partial_sum, unsigned long mask, int r[n])
{
    /* foi descoberta a soma */
    if (partial_sum == desired_sum)
    {
        /* guardar a máscara no array b */
        for (int i = 0; i < n; i++)
            r[i] = ((mask & (1 << (n-i-1))) == 0) ? 0 : 1;
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
 *
 */
int brute_force_clever(int n, integer_t p[n], integer_t desired_sum, int current_index, integer_t partial_sum, unsigned long mask, int r[n])
{
    /* foi descoberta a soma */
    if (partial_sum == desired_sum)
    {
        /* guardar a máscara no array b */
        for (int i = 0; i < n; i++)
            r[i] = ((mask & (1 << (n-i-1))) == 0) ? 0 : 1;
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

/**
 * @brief Determina a combinação dos valores de p cujo somatório é desired_sum, pelo método de Horowitz e Sahni.
 *
 * @param n                 Tamanho de p
 * @param p                 Conjunto com os valores a somar
 * @param desired_sum       Soma desejada
 * @param r                 Array que guarda a solução
 *
 * @return                  1 se foi encontrada uma solução, 0 caso contrário
 *
 */
int horowitz_and_sahni(int n, integer_t p[n], integer_t desired_sum, int r[n])
{
    int size_a = n / 2;
    int size_b = n - size_a;

    int size_a_sum = 1 << size_a;
    int size_b_sum = 1 << size_b;

    hs_data_t a[size_a_sum], b[size_b_sum];

    /* somas para o array a */
    for (int i = 0; i < size_a_sum; i++)
    {
        a[i].mask = i;
        a[i].sum = 0;
        for (int j = 0; j < n; j++)
            if ((i & (1 << (size_a-j-1))) != 0)
                a[i].sum += p[j];
    }

    /* somas para o array b */
    for (int i = 0; i < size_b_sum; i++)
    {
        b[i].mask = i;
        b[i].sum = 0;
        for (int j = 0; j < n; j++)
            if ((i & (1 << (size_b-j-1))) != 0)
                b[i].sum += p[j+size_a-1];
    }

    /* sort do array a */
    qsort(a, size_a, sizeof(hs_data_t), hs_data_cmpfunc);

    /* sort do array b */

    for (int i = 0; i < size_a_sum; i++)
        printf("%lld ",a[i].sum);
    printf("\n");
    for (int i = 0; i < size_b_sum; i++)
        printf("%lld ",b[i].sum);
    printf("\n");

    return 0;
}

/* ------------------------------------- Programa Principal ------------------------------------- */

int main(void)
{
    fprintf(stderr, "Program configuration:\n");
    fprintf(stderr, "  min_n ....... %d\n", min_n);
    fprintf(stderr, "  max_n ....... %d\n", max_n);
    fprintf(stderr, "  n_sums ...... %d\n", n_sums);
    fprintf(stderr, "  n_problems .. %d\n", n_problems);
    fprintf(stderr, "  integer_t ... %d bits\n", 8 * (int)sizeof(integer_t));

    int ret;       // valor de retorno das funções
    double t1, t2; // instantes inicial e final

#if !DEBUG
    /* cabeçalho da tabela */
    printf("n\t");
    for (int i = 1; i <= 20; i++)
        printf("t%d\t", i);
    printf("\n");
#endif

    /*
     * para cada problema
     */
    for (int i = 0; i < n_problems; i++)
    {
        int n = all_subset_sum_problems[i].n;        // número de valores a somar
        integer_t *p = all_subset_sum_problems[i].p; // conjunto com os valores a somar
        int r[n];                                    // resultado

        /* ignorar problemas com n superior */
        if (n > N_LIMIT)
            continue;

            /*
             * impressão na consola do n do problema atual
             */
#if DEBUG
        printf("n=%d\n", n);
#else
        printf("%d\t", n);
#endif

        /*
         * para cada soma
         */
        for (int j = 0; j < n_sums; j++)
        {
            integer_t desired_sum = all_subset_sum_problems[i].sums[j]; // soma desejada

            /* instante inicial */
            t1 = cpu_time();

            /*
             * execução da função
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
            ret = schroeppel_and_shamir();
#endif

            /* instante final */
            t2 = cpu_time();

            /*
             * impressão na consola
             */
#if DEBUG
            /* valor de retorno */
            printf("ret=%d\t", ret);
            /* tempo de execução */
            printf("%f\t", t2 - t1);
            /* resultado */
            if (ret == 1)
            {
                for (int i = 0; i < n; i++)
                    printf("%d", r[i]);
                printf("\n");
            }
            else
            {
                printf("Não foi encontrada uma solução!\n");
            }
#else
            /* tempo de execução */
            (ret == 1) ? printf("%f\t", t2 - t1) : printf("NaN\t");
#endif
        }

        printf("\n");
    }

    return 0;
}
