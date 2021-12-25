/*
 *
 * João Fonseca, Dezembro 2021
 *
 * Implementação das min e max heaps binárias para o método de Schroeppel e Shamir do primeiro trabalho prático de AED
 *
 */

/* ----------------------------------- Inclusão de Ficheiros ------------------------------------ */

#include <stdlib.h>
#include <stdio.h>

/* ------------------------------------ Estruturas de Dados ------------------------------------- */

typedef unsigned long long integer_t;   // 64-bit unsigned integer

typedef struct
{
    integer_t sum;
    unsigned int i1;
    unsigned int i2;
} heap_data_t;

typedef struct
{
    unsigned int capacity;
    unsigned int count;
    unsigned int type; // 0 - min-heap, 1 - max-heap
    heap_data_t *data;
} heap_t;

/* ------------------------------------- Funções de Suporte ------------------------------------- */

/**
 * @brief Retorna o índice do pai.
 * 
 * @param i Índice do elemento
 * 
 * @return  Índice do pai 
 */
unsigned int parent(unsigned int i) { return (i - 1) / 2; }

/**
 * @brief Retorna o índice do filho esquerdo.
 * 
 * @param i Índice do elemento
 * 
 * @return  Índice do filho esquerdo
 */
unsigned int left_child(unsigned int i) { return 2 * i + 1; }

/**
 * @brief Retorna o índice do filho direito.
 * 
 * @param i Índice do elemento
 * 
 * @return  Índice do filho direito 
 */
unsigned int right_child(unsigned int i) { return 2 * i + 2; }

/**
 * @brief Troca dois elementos da heap.
 * 
 * @param h     Ponteiro para a heap
 * @param i1    Índice do primeiro elemento
 * @param i2    Índice do segundo elemento
 */
void swap(heap_t *h, unsigned int i1, unsigned int i2)
{
    heap_data_t tmp = h->data[i1];
    h->data[i1] = h->data[i2];
    h->data[i2] = tmp;
}

/**
 * @brief Mantém as propriedades da heap, de baixo para cima.
 * 
 * @param h Ponteiro para a heap
 * @param i Índice do elemento
 */
void sift_up(heap_t *h, unsigned int i)
{
    unsigned int p = (i == 0) ? 0 : parent(i); // índice do pai

    /**
     * h é uma min-heap e a soma do pai é superior
     * OU
     * h é uma max-heap e a soma do pai é inferior
     */
    if ((h->type == 0 && h->data[p].sum > h->data[i].sum) || (h->type == 1 && h->data[p].sum < h->data[i].sum))
    {
        swap(h, i, p); // trocar os dois elementos
        sift_up(h, p); // fazer o mesmo para o pai
    }
}

/**
 * @brief Mantém as propriedades da heap, de cima para baixo.
 * 
 * @param h Ponteiro para a heap
 * @param i Índice do elemento
 */
void sift_down(heap_t *h, unsigned int i)
{
    unsigned int lc = left_child(i);  // índice do filho esquerdo
    unsigned int rc = right_child(i); // índice do filho direito

    lc = (lc >= h->count) ? 0 : lc; // lc = 0 se não existir filho esquerdo
    rc = (rc >= h->count) ? 0 : rc; // rc = 0 se não existir filho direito

    /* i não tem filhos */
    if (lc == rc)
        return;

    /**
     * i apenas tem lc                                      //
     * OU                                                   //
     * h é uma min-heap e a soma do lc é inferior à do rc   // chosen_c = lc, caso contrário chosen_c = rc
     * OU                                                   // 
     * h é uma max-heap e a soma do lc é superior à do rc   //
     */
    unsigned int chosen_c = ((rc < lc) || (h->type == 0 && h->data[lc].sum < h->data[rc].sum) || (h->type == 1 && h->data[lc].sum > h->data[rc].sum)) ? lc : rc;

    /**
     * h é uma min-heap e a soma do chosen_c é inferior à do i
     * OU
     * h é uma max-heap e a soma do chosen_c é superior à do i 
     */
    if ((h->type == 0 && h->data[chosen_c].sum < h->data[i].sum) || (h->type == 1 && h->data[chosen_c].sum > h->data[i].sum))
    {
        swap(h, i, chosen_c);
        sift_down(h, chosen_c);
    }
}

/* ------------------------------------- Funções Principais ------------------------------------- */

/**
 * @brief Cria uma heap.
 *
 * @param capacity  Número máximo de elementos na heap
 * @param type      Tipo de heap (0 - min-heap, 1 - max-heap)
 *
 * @return          Ponteiro para a heap se for criada com sucesso, ponteiro para NULL caso contrário
 */
heap_t *create(unsigned int capacity, unsigned int type)
{
    heap_t *heap = (heap_t *)malloc(sizeof(heap_t));

    if (heap == NULL)
        return NULL;

    heap->capacity = capacity;
    heap->count = 0;
    heap->type = type;
    heap->data = (heap_data_t *)malloc(capacity * sizeof(heap_data_t));

    if (heap->data == NULL)
        return NULL;

    return heap;
}

/**
 * @brief Destroi uma heap.
 *
 * @param h Ponteiro para a heap
 */
void destroy(heap_t *h)
{
    free(h->data);
    free(h);
}

/**
 * @brief Insere um elemento numa heap.
 *
 * @param h     Ponteiro para a heap
 * @param data  Elemento
 */
void push(heap_t *h, heap_data_t data)
{
    /* a heap tem espaço livre */
    if (h->count < h->capacity)
    {
        h->data[h->count] = data; // insere o elemento no fim da heap
        sift_up(h, h->count);     // garante as propriedades da heap
        h->count++;               // incrementa o número de elementos na heap
    }
}

/**
 * @brief Retorna o elemento no topo da heap.
 *
 * @param h Ponteiro para a heap
 *
 * @return  Elemento
 */
heap_data_t peek(heap_t *h)
{
    return h->data[0];
}

/**
 * @brief Remove e retorna o elemento no topo da heap.
 *
 * @param h Ponteiro para a heap
 * @return  Elemento
 */
heap_data_t pop(heap_t *h)
{
    heap_data_t data = peek(h);

    swap(h, 0, h->count - 1); // substitui o primeiro elemento pelo último
    h->count--;               // decrementa o número de elementos na heap
    sift_down(h, 0);          // garante as propriedades da heap

    return data;
}