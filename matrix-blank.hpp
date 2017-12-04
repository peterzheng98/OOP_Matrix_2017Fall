#ifndef SJTU_MATRIX_HPP
#define SJTU_MATRIX_HPP

#include <cstddef>
#include <initializer_list>
#include <utility>
#include <iterator>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <cassert>

using std::size_t;

namespace sjtu {

    template<class T>
    class Matrix {
    private:
        // your private member variables here.
        std::vector<std::vector<T> > _matrix;
        size_t _row, _col;
    public:
        Matrix() = default;

        Matrix(std::vector<std::vector<T> > __matrix, unsigned int __row, unsigned int __col){
            _matrix = __matrix, _row = __row, _col = __col;
        }

        Matrix(std::vector<std::vector<T> > __matrix, size_t __row, size_t __col){
            _matrix = __matrix, _row = __row, _col = __col;
        }
        Matrix matrixData(){
            return _matrix;
        }
        Matrix(size_t n, size_t m, T _init = T()) {
            std::invalid_argument e("length_error");
            //assert(n <= 0 || m <= 0);
            if (n <= 0 || m <= 0) throw e;
            std::vector<std::vector<T> >().swap(_matrix);
            _matrix.resize(n);
            for (int i = 0; i < n; i++) {
                std::vector<T> p;
                for (int j = 0; j < m; j++) {
                    p.push_back(_init);
                }
                _matrix.push_back(p);
            }
            _row = n;
            _col = m;
        }

        explicit Matrix(std::pair<size_t, size_t> sz, T _init = T()) {
            std::invalid_argument e("length_error");
            //assert(sz.first <= 0 || sz.second <= 0);
            if (sz.first <= 0 || sz.second <= 0) throw e;
            std::vector<std::vector<T> >().swap(_matrix);
            _matrix.resize(sz.first);
            for (int i = 0; i < sz.first; i++) {
                std::vector<T> p;
                for (int j = 0; j < sz.second; j++) {
                    p.push_back(_init);
                }
                _matrix.push_back(p);
            }
        }

        Matrix(const Matrix &o) {
            _row = o.rowLength(), _col = o.columnLength();
            _matrix = o.matrixData();
        }

        template<class U>
        Matrix(const Matrix<U> &o) {
            _row = o.rowLength(), _col = o.columnLength();
            std::vector<std::vector<U> > __matrix;
            for (int i = 0; i < _row; i++) {
                std::vector<U> p;
                for (int j = 0; j < _col; j++) {
                    p.push_back((T) o(i,j));
                }
                __matrix.push_back(p);
            }
            _matrix = __matrix;
        }

        Matrix &operator=(const Matrix &o) {
            _row = o.rowLength(), _col = o.columnLength();
            _matrix = o.matrixData();
            return *this;
        }

        template<class U>
        Matrix &operator=(const Matrix<U> &o) {
            _row = o.rowLength(), _col = o.columnLength();
            std::vector<std::vector<T> >().swap(_matrix);
            for (int i = 0; i < _row; i++) {
                std::vector<T> p;
                for (int j = 0; j < _col; j++) {
                    p.push_back((T) o(i,j));
                }
                _matrix.push_back(p);
            }
            return *this;
        }

        Matrix(Matrix &&o) noexcept {
            _row = o.rowLength(), _col = o.columnLength();
            std::vector<std::vector<T> >().swap(_matrix);
            for (int i = 0; i < _row; i++) {
                std::vector<T> p;
                for (int j = 0; j < _col; j++) {
                    p.push_back((T) o(i,j));
                }
                _matrix.push_back(p);
            }
        }

        Matrix &operator=(Matrix &&o) noexcept {
            _row = o.rowLength(), _col = o.columnLength();
            _matrix = o.matrixData();
            return *this;
        }

        ~Matrix() {
            //std::vector<std::vector<T> >().swap(_matrix);
        }

        Matrix(std::initializer_list<std::initializer_list<T>> il) {
            std::invalid_argument e("length_error");
            //assert(il.size() <= 0);
            if (il.size() <= 0) throw e;
            std::vector<std::vector<T> >().swap(_matrix);
            for (auto it : il) {
                std::vector<T> p;
                for (auto _it : it) {
                    p.push_back(_it);
                }
                _matrix.push_back(p);
            }
        }

    public:
        size_t rowLength() const {
            return _row;
        }

        size_t columnLength() const {
            return _col;
        }

        void resize(size_t _n, size_t _m, T _init = T()) {
            std::invalid_argument e("length_error");
            //assert(_n <= 0 || _m <= 0);
            if (_n <= 0 || _m <= 0) throw e;
            std::vector<std::vector<T> > __matrix;
            for (int i = 0; i < _n; i++) {
                std::vector<T> p;
                for (int j = 0; j < _m; j++) {
                    if (i >= _row) p.push_back(_init);
                    else if (j >= _col) p.push_back(_init);
                    else p.push_back(_matrix.at(i).at(j));
                }
                __matrix.push_back(p);
            }
            _matrix = __matrix;
        }

        void resize(std::pair<size_t, size_t> sz, T _init = T()) {
            std::invalid_argument e("length_error");
            //assert(sz.first<=0||sz.second<=0);
            if (sz.first <= 0 || sz.second <= 0) throw e;
            std::vector<std::vector<T> > __matrix;
            for (int i = 0; i < sz.first; i++) {
                std::vector<T> p;
                for (int j = 0; j < sz.second; j++) {
                    if (i >= _row) p.push_back(_init);
                    else if (j >= _col) p.push_back(_init);
                    else p.push_back(_matrix.at(i).at(j));
                }
                __matrix.push_back(p);
            }
            _matrix = __matrix;
        }

        std::pair<size_t, size_t> size() const {
            return std::make_pair(__row, __col);
        };

        void clear() {
            std::vector<std::vector<T> >().swap(_matrix);
        }

    public:
        const T &operator()(size_t i, size_t j) const {
            return _matrix.at(i).at(j);
        }

        T &operator()(size_t i, size_t j) {
            return _matrix.at(i).at(j);
        }

        Matrix<T> row(size_t _i) const {
            std::vector<std::vector<T> > __matrix;
            std::vector<T> __row;
            for (int i = 0; i < _col; i++) {
                __row.push_back(_matrix[_i][i]);
            }
            __matrix.push_back(__row);
            Matrix<T> ___matrix(__matrix, 1, _col);
            /*___matrix._col = _col;
            ___matrix._row = 1;
            ___matrix._matrix = __matrix;*/
            return ___matrix;
        }

        Matrix<T> column(size_t _i) const {
            std::vector<std::vector<T> > __matrix;
            for (int i = 0; i < _row; i++) {
                std::vector<T> __row;
                __row.push_back(_matrix[i][_i]);
                __matrix.push_back(__row);
            }
            Matrix<T> ___matrix(__matrix, _row, 1);
            /*___matrix._col = 1;
            ___matrix._row = _row;
            ___matrix._matrix = __matrix;*/
            return ___matrix;
        }


    public:
        template<class U>
        bool operator==(const Matrix<U> &o) const {
            if(o.columnLength() != _col || o.rowLength() != _row) return false;
            for(int i=0;i<_col;i++){
                for(int j=0;j<_row;j++){
                    if(o(i,j) != _matrix[i][j]) return false;
                }
            }
            return true;
        }

        template<class U>
        bool operator!=(const Matrix<U> &o) const {
            return !(o==(*this));
        }

        Matrix operator-() const {
            for(int i=0;i<_col;i++){
                for(int j=0;j<_row;j++){
                    _matrix[i][j] = -_matrix[i][j];
                }
            }
        }

        template<class U>
        Matrix &operator+=(const Matrix<U> &o) {
            //assert(_row == o._row && _col == o._col);
            if (_row == o.rowLength() && _col == o.columnLength()) {
                
                std::vector< std::vector<decltype(T() + U())> > __matrix; 
                for (int i = 0; i < _row; i++) {
                    std::vector<decltype(T() + U())> p;
                    for (int j = 0; j < _col; j++) {
                        p.push_back(_matrix[i][j] + o(i,j));
                    }
                    __matrix.push_back(p);
                }
                Matrix<decltype(T() + U())> fuck(__matrix, _row, _col);
            } else {
                std::invalid_argument e("length_error");
                throw e;
            }
            return fuck;
        }

        template<class U>
        Matrix &operator-=(const Matrix<U> &o) {
            //assert(_row == o._row && _col == o._col);
            if (_row == o.rowLength() && _col == o.columnLength()) {
                
                std::vector< std::vector<decltype(T() - U())> > __matrix; 
                for (int i = 0; i < _row; i++) {
                    std::vector<decltype(T() - U())> p;
                    for (int j = 0; j < _col; j++) {
                        p.push_back(_matrix[i][j] - o(i,j));
                    }
                    __matrix.push_back(p);
                }
                Matrix<decltype(T() - U())> fuck(__matrix, _row, _col);
            } else {
                std::invalid_argument e("length_error");
                throw e;
            }
            return fuck;
        }

        template<class U>
        Matrix &operator*=(const U &x) {
            std::invalid_argument e("length_error");
            if(_col != x.rowLength() || x.columnLength()<=0 || x.rowLength()<=0) {
                throw e;
            } else {
            std::vector<std::vector<decltype(T() * U())> > __matrix;
            __matrix.resize(_row);
            for(int i=0;i<_row;i++){
                std::vector<T> d;
                for(int j=0;j<x.columnLength();j++){
                    decltype(T() * U()) sum = 0;
                    for(int k=0;k<_col;k++){
                        sum+=_matrix[i][k] * x(k,j);
                    }
                    d.push_back(sum);
                }
                __matrix.push_back(d);
            }
            Matrix<decltype(T()*U())> fuck_oop(__matrix, _row, x.columnLength());
            return fuck_oop;
            }
        }

        Matrix tran() const {
            std::invalid_argument e("length_error");
            if(_col <= 0 || _row <= 0)
                throw e;
            else {
                std::vector< std::vector<T> > __matrix;
                for(int i=0;i<_row;i++){
                    std::vector<T> p;
                    for(int j=0;j<_col;j++){
                        p.push_back(_matrix[j][i]);
                    }
                    __matrix.push_back(p);
                }
                Matrix<T> fuck_trans(__matrix, _col, _row);
                return fuck_trans;
            }
        }

    public: // iterator
        class iterator {
        public:
            using iterator_category = std::random_access_iterator_tag;
            using value_type        = T;
            using pointer           = T *;
            using reference         = T &;
            using size_type         = size_t;
            using difference_type   = std::ptrdiff_t;

            iterator() = default;

            iterator(const iterator &) = default;

            iterator &operator=(const iterator &) = default;

        private:


        public:
            difference_type operator-(const iterator &o) {

            }

            iterator &operator+=(difference_type offset) {

            }

            iterator operator+(difference_type offset) const {

            }

            iterator &operator-=(difference_type offset) {

            }

            iterator operator-(difference_type offset) const {

            }

            iterator &operator++() {

            }

            iterator operator++(int) {

            }

            iterator &operator--() {

            }

            iterator operator--(int) {

            }

            reference operator*() const {

            }

            pointer operator->() const {

            }

            bool operator==(const iterator &o) const {

            }

            bool operator!=(const iterator &o) const {

            }
        };

        iterator begin() {

        }

        iterator end() {

        }

        std::pair<iterator, iterator> subMatrix(std::pair<size_t, size_t> l, std::pair<size_t, size_t> r) {

        }
    };

}

//
namespace sjtu {
    template<class T, class U>
    auto operator*(const Matrix<T> &mat, const U &x) {
        std::vector< std::vector<decltype(T() * U())> > __matrix;
        for(int i=0;i<mat.rowLength();i++){
            std::vector<decltype(T() * U())> p;
            for(int j=0;j<mat.columnLength();j++)
                p.push_back(mat(i,j) * x);
            __matrix.push_back(p);
        }
        Matrix<decltype(T() * U())> ans(__matrix, mat.rowLength(), mat.columnLength());
        return ans;
    }

    template<class T, class U>
    auto operator*(const U &x, const Matrix<T> &mat) {
        return mat*U;
    }

    template<class U, class V>
    auto operator*(const Matrix<U> &a, const Matrix<V> &b) {
        std::invalid_argument e("length_error");
            if(a.columnLength() != b.rowLength() || a.columnLength()<=0 || a.rowLength()<=0 || b.rowLength() <=0 || b.columnLength()) {
                throw e;
            } else {
            std::vector<std::vector<decltype(U() * V())> > __matrix;
            for(int i=0;i<a.rowLength();i++){
                std::vector<T> d;
                for(int j=0;j<b.columnLength();j++){
                    decltype(T() * U()) sum = 0;
                    for(int k=0;k<a.columnLength();k++){
                        sum+=a(i,k) * b(k,j);
                    }
                    d.push_back(sum);
                }
                __matrix.push_back(d);
            }
            Matrix<decltype(U()*V())> fuck_oop(__matrix, a.rowLength(), b.columnLength());
            return fuck_oop;
    }

    template<class U, class V>
    auto operator+(const Matrix<U> &a, const Matrix<V> &b) {
        if(a.columnLength() != b.columnLength() || a.rowLength() != b.rowLength() || a.rowLength()<=0 || a.columnLength()<=0){
            std::invalid_argument e("length_error");
            throw e;
        } else {
            std::vector<std::vector<decltype(U() + V())> > __matrix;
            for(int i=0;i<a.rowLength;i++){
                std::vector<decltype(U() + V())> p;
                for(int j=0;j<a.columnLength();j++){
                    p.push_back(a(i,j) + b(i,j));
                }
                __matrix.push_back(p);
            }
            Matrix<decltype(U() + V())> fuckit(__matrix, a.rowLength(), a.columnLength());
            return fuckit;
        }
    }

    template<class U, class V>
    auto operator-(const Matrix<U> &a, const Matrix<V> &b) {
        if(a.columnLength() != b.columnLength() || a.rowLength() != b.rowLength() || a.rowLength()<=0 || a.columnLength()<=0){
            std::invalid_argument e("length_error");
            throw e;
        } else {
            std::vector<std::vector<decltype(U() - V())> > __matrix;
            for(int i=0;i<a.rowLength;i++){
                std::vector<decltype(U() - V())> p;
                for(int j=0;j<a.columnLength();j++){
                    p.push_back(a(i,j) - b(i,j));
                }
                __matrix.push_back(p);
            }
            Matrix<decltype(U() + V())> fuckit(__matrix, a.rowLength(), a.columnLength());
            return fuckit;
    }

}

#endif //SJTU_MATRIX_HPP

