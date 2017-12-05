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
        template<class V, class U>
        friend auto operator*(const Matrix<V> &, const U &);

        template<class U, class V>
        friend auto operator*(const U &, const Matrix<V> &);

        template<class U, class V>
        friend auto operator*(const Matrix<U> &, const Matrix<V> &);

        template<class U, class V>
        friend auto operator+(const Matrix<U> &, const Matrix<V> &);

        template<class U, class V>
        friend auto operator-(const Matrix<U> &, const Matrix<V> &);

    private:
        // your private member variables here.
        std::vector<std::vector<T> > _matrix;
        size_t _row, _col;
    public:
        Matrix() = default;

        Matrix(std::vector<std::vector<T> > __matrix, unsigned int __row, unsigned int __col) {
            _matrix = __matrix, _row = __row, _col = __col;
        }

        Matrix(std::vector<std::vector<T> > __matrix, size_t __row, size_t __col) {
            _matrix = __matrix, _row = __row, _col = __col;
        }

        Matrix(size_t n, size_t m, T _init = T()) {
            std::invalid_argument e("length_error");
            //assert(n <= 0 || m <= 0);
            if (n <= 0 || m <= 0) throw e;
            std::vector<std::vector<T> >().swap(_matrix);
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
            std::vector<std::vector<T> >().swap(_matrix);
            for (int i = 0; i < _row; i++) {
                std::vector<T> p;
                for (int j = 0; j < _col; j++)
                    p.push_back(o(i, j));
                _matrix.push_back(p);
            }
        }

        template<class U>
        Matrix(const Matrix<U> &o) {
            _row = o.rowLength();
            _col = o.columnLength();
            std::vector<std::vector<T> >().swap(_matrix);
            for (int i = 0; i < _row; i++) {
                std::vector<T> p;
                for (int j = 0; j < _col; j++) {
                    p.push_back(T(o(i, j)));
                }
                _matrix.push_back(p);
            }
        }

        Matrix &operator=(const Matrix &o) {
            if (&o != this) {
                this->_col = o.columnLength();
                this->_row = o.rowLength();
                this->_matrix.resize(_row);
                for (int i = 0; i < this->_row; i++)
                    this->_matrix[i].resize(_col);

                for (int i = 0; i < o.rowLength(); i++) {
                    for (int j = 0; j < o.columnLength(); j++)
                        this->_matrix[i][j] = o(i, j);
                }
            }
            return *this;
        }

        template<class U>
        Matrix &operator=(const Matrix<U> &o) {
            _row = o.rowLength(), _col = o.columnLength();
//            std::vector<std::vector<T> > newp;
            std::vector<std::vector<T> >().swap(_matrix);
            for (int i = 0; i < _row; i++) {
                std::vector<T> p;
                for (int j = 0; j < _col; j++) {
                    p.push_back(T(o(i, j)));
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
                    p.push_back((T) o(i, j));
                }
                _matrix.push_back(p);
            }
        }

        Matrix &operator=(Matrix &&o) noexcept {
            _row = o.rowLength(), _col = o.columnLength();
            std::vector<std::vector<T> >().swap(_matrix);
            for (int i = 0; i < _row; i++) {
                std::vector<T> p;
                for (int j = 0; j < _col; j++)
                    p.push_back(o(i, j));
                _matrix.push_back(p);
            }
            return *this;
        }

        ~Matrix() {
            //std::vector<std::vector<T> >().swap(_matrix);
        }

        Matrix(std::initializer_list<std::initializer_list<T>> il) {
            std::invalid_argument e("length_error");
            //assert(il.size() <= 0);
            size_t __col = 0, __row = 0;
            if (il.size() <= 0) throw e;
            std::vector<std::vector<T> >().swap(_matrix);
            for (auto it : il) {
                std::vector<T> p;
                __col = 0;
                for (auto _it : it) {
                    p.push_back(_it);
                    __col++;
                }
                _matrix.push_back(p);
                __row++;
            }
            _col = __col, _row = __row;
        }

    public:
        size_t rowLength() const {
            return _matrix.size();
        }

        size_t columnLength() const {
            return _matrix.size() > 0 ? _matrix[0].size() : 0;
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
            return std::make_pair(_matrix.size(), _matrix.size() > 0 ? _matrix[0].size() : 0);
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
            if (o.columnLength() != _col || o.rowLength() != _row) return false;
            for (int i = 0; i < _row; i++) {
                for (int j = 0; j < _col; j++) {
                    if (o(i, j) != _matrix[i][j]) return false;
                }
            }
            return true;
        }

        template<class U>
        bool operator!=(const Matrix<U> &o) const {
            return !(o == (*this));
        }

        Matrix operator-() const {
            std::vector<std::vector<T> > __matrix;
            for (int i = 0; i < _row; i++) {
                std::vector<T> p;
                for (int j = 0; j < _col; j++) {
                    p.push_back(-_matrix[i][j]);
                }
                __matrix.push_back(p);
            }
            Matrix<T> ans(__matrix, _row, _col);
            return ans;
        }

        template<class U>
        Matrix &operator+=(const Matrix<U> &o) {
            //assert(_row == o._row && _col == o._col);
            if (_row == o.rowLength() && _col == o.columnLength()) {

                std::vector<std::vector<decltype(T() + U())> > __matrix;
                for (int i = 0; i < _row; i++) {
                    std::vector<decltype(T() + U())> p;
                    for (int j = 0; j < _col; j++) {
                        p.push_back(_matrix[i][j] + o(i, j));
                    }
                    __matrix.push_back(p);
                }
                Matrix<decltype(T() + U())> fuck(__matrix, _row, _col);
                *this = fuck;
                return *this;
            } else {
                std::invalid_argument e("length_error");
                throw e;
            }
        }

        template<class U>
        Matrix &operator-=(const Matrix<U> &o) {
            //assert(_row == o._row && _col == o._col);
            if (_row == o.rowLength() && _col == o.columnLength()) {

                std::vector<std::vector<decltype(T() - U())> > __matrix;
                for (int i = 0; i < _row; i++) {
                    std::vector<decltype(T() - U())> p;
                    for (int j = 0; j < _col; j++) {
                        p.push_back(_matrix[i][j] - o(i, j));
                    }
                    __matrix.push_back(p);
                }
                Matrix<decltype(T() - U())> fuck(__matrix, _row, _col);
                *this = fuck;
                return *this;
            } else {
                std::invalid_argument e("length_error");
                throw e;
            }
        }

        template<class U>
        Matrix &operator*=(const U &x) {
            std::vector<std::vector<decltype(T() * U())> > __matrix;
            for (int i = 0; i < this->_row; i++) {
                std::vector<decltype(T() * U())> d;
                for (int j = 0; j < this->_col; j++) {
                    d.push_back(this->_matrix[i][j] * x);
                }
                __matrix.push_back(d);
            }
            Matrix<decltype(T() * U())> fuck_oop(__matrix, this->_row, this->_col);
            *this = fuck_oop;
            return *this;
        }


        Matrix tran() const {
            std::invalid_argument e("length_error");
            if (_col <= 0 || _row <= 0)
                throw e;
            else {
                std::vector<std::vector<T> > __matrix;
                for (int i = 0; i < _col; i++) {
                    std::vector<T> p;
                    for (int j = 0; j < _row; j++) {
                        p.push_back(_matrix[j][i]);
                    }
                    __matrix.push_back(p);
                }
                Matrix<T> fuck_trans(__matrix, __matrix.size(), __matrix[0].size());
                return fuck_trans;
//                std::vector<std::vector<T> > __matrix;
//                Matrix fuck_trans(__matrix,columnLength(),rowLength());
//                for (size_t i = 0; i < columnLength(); i ++){
//                    for (size_t j = 0; j < rowLength(); j ++){
//                        fuck_trans(i, j) = (*this)(j, i);
//                    }
//                }
//                return fuck_trans;
            }

        }

        void print() const {
            for (int i = 0; i < _row; i++) {
                for (int j = 0; j < _col; j++)
                    printf("%d ", _matrix[i][j]);
                printf("\n");
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

            iterator(Matrix *p, size_t curI) {
                _limit = false;
                _p = p;
                _curI = curI;
            }

            iterator(Matrix *p, size_t curI, bool limit, size_t next_limit = 0, size_t cur_limit = 0,
                     size_t head_limit = 0) {
                _p = p, _curI = curI;
                _limit = limit, _next_limit = next_limit;
                _cur_limit = cur_limit, _head_limit = head_limit;
            }

            iterator(const iterator &) = default;

            iterator &operator=(const iterator &) = default;

        private:
            Matrix *_p;
            std::pair<size_t, size_t> cur;
            size_t _curI;
            bool _limit;
            size_t _next_limit, _cur_limit, _head_limit;
        public:
            difference_type operator-(const iterator &o) {
//                std::cout << "o._curI - _curI = " << o._curI - _curI << std::endl;
//                std::cout << "o._curI = " << o._curI << std::endl;
//                std::cout << "_curI = " << _curI << std::endl;
                int p1 = o._curI, p2 = _curI;
                return p2 - p1;
            }

            iterator &operator+=(difference_type offset) {
                iterator ite(_p, _curI + offset);
                *this = ite;
                return *this;
            }

            iterator operator+(difference_type offset) const {
                iterator ite(_p, _curI + offset);
                return ite;
            }

            iterator &operator-=(difference_type offset) {
                iterator ite(_p, _curI - offset);
                *this = ite;
                return *this;
            }

            iterator operator-(difference_type offset) const {
                iterator ite(_p, _curI - offset);
                return ite;
            }

            iterator &operator++() {
//                std::cout << _limit << std::endl;
                if (!_limit) {
                    *this += 1;
                    return *this;
                } else {
//                    std::cout << "RUN IN SPEC" << std::endl;
                    size_t __curI = _curI + 1;
                    if (__curI % _head_limit > _next_limit)
                        __curI = (_curI / _head_limit + 1) * _head_limit + _cur_limit;
                    if (__curI % _head_limit < _cur_limit)
                        __curI = (_curI / _head_limit + 1) * _head_limit + _cur_limit;
                    iterator ite(_p, __curI, _limit, _next_limit, _cur_limit, _head_limit);
                    *this = ite;
                    return *this;
                }

            }

            iterator operator++(int) {
                iterator cur(_p, _curI);
                *this += 1;
                return cur;
            }

            iterator &operator--() {
                *this -= 1;
                return *this;
            }

            iterator operator--(int) {
                iterator cur(_p, _curI);
                *this -= 1;
                return cur;
            }

            reference operator*() const {
                size_t R = _curI / (_p->columnLength());
                size_t C = _curI % (_p->columnLength());
                return _p->_matrix[R][C];
            }

            pointer operator->() const {
                size_t R = _curI / (_p->columnLength());
                size_t C = _curI % (_p->columnLength());
                return &(_p->_matrix[R][C]);
            }

            bool operator==(const iterator &o) const {
                return o._curI == _curI && o._p == _p;
            }

            bool operator!=(const iterator &o) const {
                return !(o._curI == _curI && o._p == _p);
            }
        };

        iterator begin() {
//            return &_matrix[start.first][start.second];
            return iterator(this, 0);
        }

        iterator end() {
//            return &(_matrix[finish.first][finish.second]+1);
            return iterator(this, _col * _row);
        }

        std::pair<iterator, iterator> subMatrix(std::pair<size_t, size_t> l, std::pair<size_t, size_t> r) {
            return std::make_pair(iterator(this, (l.first * this->columnLength() + l.second), true, r.second, l.second,
                                           this->columnLength()),
                                  iterator(this, (r.first * this->columnLength() + r.second), true, r.second, l.second,
                                           this->columnLength()));
        }
    };

}

//
namespace sjtu {
    template<class V, class U>
    auto operator*(const Matrix<V> &mat, const U &x) {
        std::vector<std::vector<decltype(V() * U())> > __matrix;
        for (int i = 0; i < mat.rowLength(); i++) {
            std::vector<decltype(U() * V())> p;
            for (int j = 0; j < mat.columnLength(); j++)
                p.push_back(mat._matrix[i][j] * x);
            __matrix.push_back(p);
        }
        Matrix<decltype(V() * U())> ans(__matrix, mat.rowLength(), mat.columnLength());
        return ans;
        /*Matrix<decltype(V() * U())> ans(mat.rowLength(), mat.columnLength());
//        std::vector<std::vector<decltype(V() * U())> > __matrix;
        for (int i = 0; i < mat.rowLength(); i++) {
//            std::vector<decltype(V() * U())> p;
            for (int j = 0; j < mat.columnLength(); j++)
                ans._matrix[i][j] = mat._matrix[i][j] * x;
//            __matrix.push_back(p);
        }
//        Matrix<decltype(V() * U())> ans(__matrix, mat.rowLength(), mat.columnLength());
        return ans;*/
    }

    template<class U, class V>
    auto operator*(const U &x, const Matrix<V> &mat) {
        std::vector<std::vector<decltype(U() * V())> > __matrix;
        for (int i = 0; i < mat.rowLength(); i++) {
            std::vector<decltype(U() * V())> p;
            for (int j = 0; j < mat.columnLength(); j++)
                p.push_back(mat._matrix[i][j] * x);
            __matrix.push_back(p);
        }
        Matrix<decltype(U() * V())> ans(__matrix, mat.rowLength(), mat.columnLength());
        return ans;


    }

    template<class U, class V>
    auto operator*(const Matrix<U> &a, const Matrix<V> &b) {
        std::invalid_argument e("length_error");
        if (a.columnLength() != b.rowLength() || a.columnLength() < 0 || a.rowLength() < 0 || b.rowLength() < 0 ||
            b.columnLength() < 0) {
            throw e;
        } else {
            std::vector<std::vector<decltype(U() * V())> > __matrix;
            for (int i = 0; i < a.rowLength(); i++) {
                std::vector<decltype(U() * V())> d;
                for (int j = 0; j < b.columnLength(); j++) {
                    decltype(V() * U()) sum = 0;
                    for (int k = 0; k < a.columnLength(); k++) {
                        sum += a._matrix[i][k] * b._matrix[k][j];
                    }
                    d.push_back(sum);
                }
                __matrix.push_back(d);
            }
            Matrix<decltype(U() * V())> fuck_oop(__matrix, a.rowLength(), b.columnLength());
            return fuck_oop;
        }
    }

    template<class U, class V>
    auto operator+(const Matrix<U> &a, const Matrix<V> &b) {
        if (a.columnLength() != b.columnLength() || a.rowLength() != b.rowLength() || a.rowLength() <= 0 ||
            a.columnLength() <= 0) {
            std::invalid_argument e("length_error");
            throw e;
        } else {
            std::vector<std::vector<decltype(U() + V())> > __matrix;
            Matrix<decltype(U() + V())> fuckit(a.rowLength(), a.columnLength());
            for (int i = 0; i < a.rowLength(); i++) {
//                std::vector<decltype(U() + V())> p;
                for (int j = 0; j < a.columnLength(); j++) {
                    fuckit._matrix[i][j] = a._matrix[i][j] + b._matrix[i][j];
                }
//                __matrix.push_back(p);
            }

            return fuckit;
        }
    }

    template<class U, class V>
    auto operator-(const Matrix<U> &a, const Matrix<V> &b) {
        if (a.columnLength() != b.columnLength() || a.rowLength() != b.rowLength() || a.rowLength() <= 0 ||
            a.columnLength() <= 0) {
            std::invalid_argument e("length_error");
            throw e;
        } else {
            std::vector<std::vector<decltype(U() + V())> > __matrix;
            Matrix<decltype(U() + V())> fuckit(a.rowLength(), a.columnLength());
            for (int i = 0; i < a.rowLength(); i++) {
//                std::vector<decltype(U() + V())> p;
                for (int j = 0; j < a.columnLength(); j++) {
                    fuckit._matrix[i][j] = a._matrix[i][j] - b._matrix[i][j];
                }
//                __matrix.push_back(p);
            }

            return fuckit;
        }
    }
}

#endif //SJTU_MATRIX_HPP

