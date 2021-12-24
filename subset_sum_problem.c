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
#define N_LIMIT 12
#endif

/*
 * 0 - brute force não recursiva
 * 1 - brute force recursiva
 * 2 - brute force clever
 * 3 - horowitz and sahni
 * 4 - schroeppel and shamir
 */
#ifndef FUNC
#define FUNC 4
#endif

/* ------------------------------------ Estruturas de Dados ------------------------------------- */

typedef struct
{
    integer_t sum;
    unsigned long long mask;
} hs_data_t;

typedef struct
{
    integer_t sum;
    unsigned long long mask;
} ss_data_t;

typedef struct
{
    integer_t sum;
    unsigned long long mask;
    unsigned int idx1;
    unsigned int idx2;
} ss_heap_data_t;


/* ------------------------------------- Funções de Suporte ------------------------------------- */

int hs_data_cmpfunc(const void *d1, const void *d2)
{
    hs_data_t *data1 = (hs_data_t *)d1;
    hs_data_t *data2 = (hs_data_t *)d2;
    return (data1->sum > data2->sum) - (data1->sum < data2->sum);
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
    unsigned long long mask; // máscara
    integer_t test_sum;      // soma de teste

    /* para cada combinação */
    for (mask = 0; mask < (1 << n); mask++)
    {
        /* determinar a soma dos valores de p */
        test_sum = 0;
        for (int i = 0; i < n; i++)
            if ((mask & (1 << (n - i - 1))) != 0)
                test_sum += p[i];

        /* foi descoberta a soma */
        if (test_sum == desired_sum)
        {
            /* guardar a máscara no array r */
            for (int i = 0; i < n; i++)
                r[i] = ((mask & (1 << (n - i - 1))) == 0) ? 0 : 1;
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
int brute_force_recursive(int n, integer_t p[n], integer_t desired_sum, int current_index, integer_t partial_sum, unsigned long long mask, int r[n])
{
    /* foi descoberta a soma */
    if (partial_sum == desired_sum)
    {
        /* guardar a máscara no array r */
        for (int i = 0; i < n; i++)
            r[i] = ((mask & (1 << (n - i - 1))) == 0) ? 0 : 1;
        return 1;
    }

    /* valor de retorno é zero se ainda não tiver sido encontrada uma solução */
    int ret = 0;
    if (current_index < n)
    {
        /* p[current_index] não é usado na soma */
        ret |= brute_force_recursive(n, p, desired_sum, current_index + 1, partial_sum, mask, r);
        /* p[current_index] é usado na soma */
        ret |= brute_force_recursive(n, p, desired_sum, current_index + 1, partial_sum + p[current_index], mask | (1 << (n - current_index - 1)), r);
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
int brute_force_clever(int n, integer_t p[n], integer_t desired_sum, int current_index, integer_t partial_sum, unsigned long long mask, int r[n])
{
    /* foi descoberta a soma */
    if (partial_sum == desired_sum)
    {
        /* guardar a máscara no array r */
        for (int i = 0; i < n; i++)
            r[i] = ((mask & (1 << (n - i - 1))) == 0) ? 0 : 1;
        return 1;
    }

    /* valor de retorno é zero se ainda não tiver sido encontrada uma solução */
    int ret = 0;
    if (current_index < n && partial_sum < desired_sum)
    {
        /* p[current_index] não é usado na soma */
        ret |= brute_force_clever(n, p, desired_sum, current_index + 1, partial_sum, mask, r);
        /* p[current_index] é usado na soma */
        ret |= brute_force_clever(n, p, desired_sum, current_index + 1, partial_sum + p[current_index], mask | (1 << (n - current_index - 1)), r);
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
    /* tamanho das metades de p */
    unsigned int size_p1 = n / 2;
    unsigned int size_p2 = n - size_p1;
    
    /* tamanho de a e b */
    unsigned long long size_a = 1 << size_p1; // 2^size_p1
    unsigned long long size_b = 1 << size_p2; // 2^size_p2

    /* alocação de a e b em memória */
    hs_data_t *a = (hs_data_t *)malloc(size_a * sizeof(hs_data_t));
    hs_data_t *b = (hs_data_t *)malloc(size_b * sizeof(hs_data_t));

    /* gerar elementos de a */
    for (unsigned long long i = 0; i < size_a; i++)
    {
        /* máscara */
        a[i].mask = i;
        /* soma */
        a[i].sum = 0;
        for (unsigned int j = 0; j < size_p1; j++)
            if (a[i].mask & (1 << (size_p1 - j - 1)))
                a[i].sum += p[j];
    }

    /* gerar elementos de b */
    for (unsigned long long i = 0; i < size_b; i++)
    {
        /* máscara */
        b[i].mask = i;
        /* soma */
        b[i].sum = 0;
        for (unsigned int j = 0; j < size_p2; j++)
            if (b[i].mask & (1 << (size_p2 - j - 1)))
                b[i].sum += p[j + size_p1];
    }

    /* ordenação (quicksort) dos elementos de a e b */
    qsort(a, size_a, sizeof(hs_data_t), hs_data_cmpfunc);
    qsort(b, size_b, sizeof(hs_data_t), hs_data_cmpfunc);

    unsigned int ret = 0;            // valor de retorno é zero se não for encontrada nenhuma solução */
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
                r[k] = (a[i].mask & (1 << (size_p1 - k - 1))) == 0 ? 0 : 1;

            /* guardar a máscara de b[j] em r */
            for (k = 0; k < size_p2; k++)
                r[k + size_p1] = (b[j].mask & (1 << (size_p2 - k - 1))) == 0 ? 0 : 1;

            ret = 1;
            break;
        }
    }

    /* desalocação de a e b em memória */
    free(a);
    free(b);

    return ret;
}

/**
 * @brief Determina a combinação dos valores de p cujo somatório é desired_sum, pelo método de Schroeppel e Shamir.
 *
 * @param n                 Tamanho de p
 * @param p                 Conjunto com os valores a somar
 * @param desired_sum       Soma desejada
 * @param r                 Array que guarda a solução
 *
 * @return                  1 se foi encontrada uma solução, 0 caso contrário
 *
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
    unsigned long long size_a1 = 1 << size_p1; // 2^size_p1
    unsigned long long size_a2 = 1 << size_p2; // 2^size_p2
    unsigned long long size_b1 = 1 << size_p3; // 2^size_p3
    unsigned long long size_b2 = 1 << size_p4; // 2^size_p4

    /* alocação de a e b em memória */
    ss_data_t *a1 = (ss_data_t *)malloc(size_a1 * sizeof(ss_data_t));
    ss_data_t *a2 = (ss_data_t *)malloc(size_a2 * sizeof(ss_data_t));
    ss_data_t *b1 = (ss_data_t *)malloc(size_b1 * sizeof(ss_data_t));
    ss_data_t *b2 = (ss_data_t *)malloc(size_b2 * sizeof(ss_data_t));
    
    /* gerar elementos de a1 */
    for (unsigned long long i = 0; i < size_a1; i++)
    {
        /* máscara */
        a1[i].mask = i;
        /* soma */
        a1[i].sum = 0;
        for (unsigned int j = 0; j < size_p1; j++)
            if (a1[i].mask & (1 << (size_p1 - j - 1)))
                a1[i].sum += p[j];
    }
    /* gerar elementos de a2 */
    for (unsigned long long i = 0; i < size_a2; i++)
    {
        /* máscara */
        a2[i].mask = i;
        /* soma */
        a2[i].sum = 0;
        for (unsigned int j = 0; j < size_p2; j++)
            if (a2[i].mask & (1 << (size_p2 - j - 1)))
                a2[i].sum += p[j + size_p1];
    }
    /* gerar elementos de b1 */
    for (unsigned long long i = 0; i < size_b1; i++)
    {
        /* máscara */
        b1[i].mask = i;
        /* soma */
        b1[i].sum = 0;
        for (unsigned int j = 0; j < size_p3; j++)
            if (b1[i].mask & (1 << (size_p3 - j - 1)))
                b1[i].sum += p[j + size_p1 + size_p2];
    }
    /* gerar elementos de b2 */
    for (unsigned long long i = 0; i < size_b2; i++)
    {
        /* máscara */
        b2[i].mask = i;
        /* soma */
        b2[i].sum = 0;
        for (unsigned int j = 0; j < size_p4; j++)
            if (b2[i].mask & (1 << (size_p4 - j - 1)))
                b2[i].sum += p[j + size_p1 + size_p2 + size_p3];
    }

    /* ordenação (quicksort) dos elementos de a1, a2, b1 e b2 */
    qsort(a1, size_a1, sizeof(ss_data_t), hs_data_cmpfunc);
    qsort(b1, size_b1, sizeof(ss_data_t), hs_data_cmpfunc);
    qsort(a2, size_a2, sizeof(ss_data_t), hs_data_cmpfunc);
    qsort(b2, size_b2, sizeof(ss_data_t), hs_data_cmpfunc);

    printf("a1=[");
    for (int i = 0; i < size_a1; i++)
        printf("%llu, ",a1[i].sum);
    printf("]\n");
    printf("a2=[");
    for (int i = 0; i < size_a2; i++)
        printf("%llu, ",a2[i].sum);
    printf("]\n");
    printf("b1=[");
    for (int i = 0; i < size_b1; i++)
        printf("%llu, ",b1[i].sum);
    printf("]\n");
    printf("b2=[");
    for (int i = 0; i < size_b2; i++)
        printf("%llu, ",b2[i].sum);
    printf("]\n");

    /* desalocação de a1, a2, b1 e b2 em memória */
    free(a1);
    free(a2);
    free(b1);
    free(b2);

    return 0;
}

/*--------------------------------------------Max-Heap--------------------------------------------*/

#define MAX_SIZE 15

// returns the index of the parent node
int parent(int i) {
    return (i - 1) / 2;
}

// return the index of the left child 
int left_child(int i) {
    return 2*i + 1;
}

// return the index of the right child 
int right_child(int i) {
    return 2*i + 2;
}

void swap(int *x, int *y) {
    int temp = *x;
    *x = *y;
    *y = temp;
}

// insert the item at the appropriate position
void insert(int a[], int data, int *n) {
    if (*n >= MAX_SIZE) {
        printf("%s\n", "The heap is full. Cannot insert");
        return;
    }
    // first insert the time at the last position of the array 
    // and move it up
    a[*n] = data;
    *n = *n + 1;


    // move up until the heap property satisfies
    int i = *n - 1;
    while (i != 0 && a[parent(i)] < a[i]) {
        swap(&a[parent(i)], &a[i]);
        i = parent(i);
    }
}

// moves the item at position i of array a
// into its appropriate position
void max_heapify(int a[], int i, int n){
    // find left child node
    int left = left_child(i);

    // find right child node
    int right = right_child(i);

    // find the largest among 3 nodes
    int largest = i;

    // check if the left node is larger than the current node
    if (left <= n && a[left] > a[largest]) {
        largest = left;
    }

    // check if the right node is larger than the current node
    if (right <= n && a[right] > a[largest]) {
        largest = right;
    }

    // swap the largest node with the current node 
    // and repeat this process until the current node is larger than 
    // the right and the left node
    if (largest != i) {
        int temp = a[i];
        a[i] = a[largest];
        a[largest] = temp;
        max_heapify(a, largest, n);
    }

}

// converts an array into a heap
void build_max_heap(int a[], int n) {
    int i;
    for (i = n/2; i >= 0; i--) {
        max_heapify(a, i, n);
    } 
}

// returns the  maximum item of the heap
int get_max(int a[]) {
    return a[0];
}

// deletes the max item and return
int extract_max(int a[], int *n) {
    int max_item = a[0];

    // replace the first item with the last item
    a[0] = a[*n - 1];
    *n = *n - 1;

    // maintain the heap property by heapifying the 
    // first item
    max_heapify(a, 0, *n);
    return max_item;
}

// prints the heap
void print_heap(int a[], int n) {
    int i;
    for (i = 0; i < n; i++) {
        printf("%d\n", a[i]);
    }
    printf("\n");
}


/* int main() {
    int n = 10;
    int a[MAX_SIZE];
    a[1] = 10; a[2] = 12; a[3] = 9; a[4] = 78; a[5] = 33; a[6] = 21; a[7] = 35; a[8] = 29; a[9] = 5; a[10] = 66;
    build_max_heap(a, n);
    insert(a, 55, &n);
    insert(a, 56, &n);
    insert(a, 57, &n);
    insert(a, 58, &n);
    insert(a, 100, &n);
    print_heap(a, n);
    return 0;
} */


/* ------------------------------------- Programa Principal ------------------------------------- */

int main(void)
{
    fprintf(stderr, "Program configuration:\n");
    fprintf(stderr, "  min_n ....... %d\n", min_n);
    fprintf(stderr, "  max_n ....... %d\n", max_n);
    fprintf(stderr, "  n_sums ...... %d\n", n_sums);
    fprintf(stderr, "  n_problems .. %d\n", n_problems);
    fprintf(stderr, "  integer_t ... %d bits\n", 8 * (int)sizeof(integer_t));
    fprintf(stderr, "  function .... %d\n", FUNC);

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
            ret = schroeppel_and_shamir(n, p, desired_sum, r);
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
