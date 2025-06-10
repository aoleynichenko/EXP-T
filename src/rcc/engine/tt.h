/*
 *  EXP-T -- A Relativistic Fock-Space Multireference Coupled Cluster Program
 *  Copyright (C) 2018-2025 The EXP-T developers.
 *
 *  This file is part of EXP-T.
 *
 *  EXP-T is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  EXP-T is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with EXP-T.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  E-mail:        exp-t-program@googlegroups.com
 *  Google Groups: https://groups.google.com/d/forum/exp-t-program
 */

/*
 * Interface to the Tensor Train Library by Artem Rumyantsev
 */

#ifndef CC_TT_H_INCLUDED
#define CC_TT_H_INCLUDED

#include <complex.h>
#include <stdint.h>


// Разложение тензора в тензорный поезд
/*
   Input:
      tensor_1d: тензор в виде одномерного массива
      len_tensor: длина массива tensor_1d
      shape: размеры тензора
      len_shape: размерность тензора, тоже самое, что длина shape
      tol: точность для сингулярного разложения
   Output:
      new_tt: тензорный поезд
*/
extern void* tt_d_new(const double *tensor_1d, const uint64_t *shape, uint64_t len_shape, double tol);
extern void* tt_z_new(const double complex *tensor_1d, const uint64_t *shape, uint64_t len_shape, double tol);


// Свертка с матрицами (singles)
/*
   Input:
      tt: тензорный поезд, с которыми сворачиваем
      tmp_ts:  массив указателей на матрицы(у), с которыми сворачиваем (однократные амплитуды)
      m_ts: кол-во строк матрицы
      n_ts: кол-во столбцов матрицы
      tt_axis: массив, который показывает по каким индексам тензора (тензорного поезда) сворачиваем
      ts_axis: массив, состоящий из 0 или 1, то есть свертка или по первому индексу, или по второму
      len_axis: сколько сверток, по сути длина массива tt_axis или ts_axis
   Output:
      new_tt: тензорный поезд
*/
extern void* tt_d_contraction_tm(const void *tt, const void *tmp_ts, uint64_t m_ts, uint64_t n_ts, const uint64_t  *tt_axis, const uint64_t *ts_axis, uint64_t len_axis);
extern void* tt_z_contraction_tm(const void *tt, const void *tmp_ts, uint64_t m_ts, uint64_t n_ts, const uint64_t  *tt_axis, const uint64_t *ts_axis, uint64_t len_axis);

// Свертка с матрицами (singles)
/*
   Input:
      tt: тензорный поезд, с которыми сворачиваем
      tmp_tt: массив указателей на тензорные поезда, с которыми сворачиваем. Размер массива может быть равен или 1, или 2. Свертка вида tmp_tt[0] tt tmp_tt[1], то есть подсоединяем с двух концов
      len_tmp_tt: размер массива
      lencontr: кол-во индексов, по которым сворачиваем
      ordered: флаг упорядочены ли индексы, данный параметр на будущее, пока возможно считать только свертки типа "ijab,bakn->ijkn"
   Output:
      new_tt: тензорный поезд
*/
extern void* tt_d_contraction_tt(const void *tt, const void *tmp_tt, uint64_t len_tmp_tt, uint64_t lencontr, _Bool ordered);
extern void* tt_z_contraction_tt(const void *tt, const void *tmp_tt, uint64_t len_tmp_tt, uint64_t lencontr, _Bool ordered);

//  Возвращение в тензор (одномерный массив)
/*
   Input:
      tt: тензорный поезд
   Output:
      tensor: тензор (одномерный массив)
*/
extern double* tt_d_to_tensor_into(const void *tt);
extern double complex* tt_z_to_tensor_into(const void *tt);

// Конвертация в массив байтов
/*
   Input:
      tt: тензорный поезд,
   Output:
      массив байтов u8
*/
extern char* tt_d_as_bytes(const void *tt, uint64_t *nbytes);
extern char* tt_z_as_bytes(const void *tt, uint64_t *nbytes);

// Конвертация из массива байтов в поезда
/*
   Input:
      tt_bytes: массив байтов
      shape_tt: одномерный массив размерностей кубиков в поезде
      nbytes: количесвто байт (для f64 8*кол-во элементов, для Complex<f64> 16*кол-во элементов)
   Output:
      Тензорный поезд
*/
extern void* tt_d_bytes_to_tt(const char *tt_bytes, uint64_t nbytes, uint64_t len_shape);
extern void* tt_z_bytes_to_tt(const char *tt_bytes, uint64_t nbytes, uint64_t len_shape);

// Перемещение индексов
/*
   Input:
      tt: тензорный поезд
      ax, bx: оси, которые меняем (только рядом)
   Output:
      Тензорный поезд
*/
extern void* tt_d_swap_near_axes(const void *tt,  uint64_t ax, uint64_t bx);
extern void* tt_z_swap_near_axes(const void *tt,  uint64_t ax, uint64_t bx);

// Печатает размеры ядер
/*
   Input:
      tt: тензорный поезд
*/
extern void tt_d_print(void* tt);
extern void tt_z_print(void* tt);

// Возвращает размеры ядер как одномерный массив
/*
   Input:
      tt: тензорный поезд
*/
extern uint64_t* tt_d_shape(void* tt);
extern uint64_t* tt_z_shape(void* tt);


// Освобождение памяти
/*
   Input:
      tt: тензорный поезд
*/

extern void tt_free(void* tt);

#endif // CC_TT_H_INCLUDED
