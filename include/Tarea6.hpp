#include "Matrix.hpp"
#include <cstdlib>
#include <vector>
#include <lapacke.h>

/**
 * Util funcion para visualizar matrices
 * @tparam T
 * @param matrix
 */
template<typename T>
void printMatrix(anpi::Matrix<T> matrix) {
    for (unsigned int i = 0; i < matrix.rows(); i++) {
        std::cout << "{ ";
        for (unsigned int j = 0; j < matrix.cols(); j++) {
            if (j == matrix.cols() - 1) {
                std::cout << matrix[i][j];
            } else std::cout << matrix[i][j] << ", ";
        }
        std::cout << " }" << std::endl;
    }
}

/**
 * Funcion que genera numeros aleatorios entre un rango de -100 a 100
 * @tparam T
 * @return
 */
template<typename T>
T getRandom() {
    T temp = T(rand() % 200 - 100); // rango -100 a  100
    return temp;
}

/**
 * Genera una matriz cuadrada llena de numeros aleatorios
 * @tparam T
 * @param N
 * @return
 */
template<typename T>
anpi::Matrix<T> randomSymmetricSqr(const size_t N) {
    anpi::Matrix<T> temp = anpi::Matrix<T>(N, N);
    for (unsigned int i = 0; i < N; i++) {
        for (unsigned int j = 0; j < N; j++) {
            temp[i][j] = getRandom<T>();
        }
    }
    return temp;
}


/**
 * Calcula un vector con los eigenvalores val y los correspondientes eigenvectores E de la matriz cuadrada simetrica A.
 * @tparam T
 * @param A
 * @param val
 * @param E
 */
template<typename T>
void eig(const anpi::Matrix<T> &A, std::vector<T> &val, anpi::Matrix<T> &E) {
    lapack_int info, LDA, N;
    LDA = A.rows();
    int x = A.rows();
    double tempa[x][x];

    for (unsigned int i = 0; i < x; i++) {
        for (unsigned int j = 0; j < x; j++) {
            tempa[i][j] = A[i][j];
        }
    }

    double tempw[x][x];
    info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', LDA, *tempa, LDA, *tempw);

    for (unsigned int i = 0; i < x; i++) {
        for (unsigned int j = 0; j < x; j++) {
            E[i][j] = tempa[i][j];
        }
    }

    for (unsigned int i = 0; i < A.rows(); i++) {
        val.push_back(tempw[0][i]);
    }

}

/**
 * Funcion de ordenamiento para el metodo de Jacobi
 * @param d eigenvalores
 * @param v matriz con los eigenvectores asociados
 */
void jacobiSort(std::vector<double> &d, anpi::Matrix<double> *v = NULL) {
    unsigned int k;
    for (unsigned int i = 0; i < d.size() - 1; i++) {
        double p = d[k = i];

        for (unsigned int j = i; j < d.size(); j++)
            if (d[j] >= p) p = d[k = j];

        if (k != i) {
            d[k] = d[i];
            d[i] = p;
            if (v != NULL)

                for (unsigned int j = 0; j < d.size(); j++) {
                    p = (*v)[j][i];
                    (*v)[j][i] = (*v)[j][k];
                    (*v)[j][k] = p;
                }
        }
    }
}

/**
 * Estructura general del Jacobi
 * Recibe la matriz del eigensistema
 */
struct JacobiStruct {

    const int n;
    anpi::Matrix<double> a, v;
    std::vector<double> d;
    int nrot;
    const double EPS;

    JacobiStruct(anpi::Matrix<double> &aa) : n(aa.rows()), a(aa), v(n, n), d(n), nrot(0),
                                  EPS(std::numeric_limits<double>::epsilon()) {
        unsigned int i, j, ip, iq;
        double tresh, theta, tau, t, sm, s, h, g, c;
        std::vector<double> b(n), z(n);
        for (ip = 0; ip < n; ip++) {
            for (iq = 0; iq < n; iq++) v[ip][iq] = 0.0;
            v[ip][ip] = 1.0;
        }

        // acumula los terminos de la ecuacion en un vector
        for (ip = 0; ip < n; ip++) {
            b[ip] = d[ip] = a[ip][ip];
            z[ip] = 0.0;
        }
        for (i = 1; i <= 50; i++) {
            sm = 0.0;
            for (ip = 0; ip < n - 1; ip++) {
                for (iq = ip + 1; iq < n; iq++)
                    sm += abs(a[ip][iq]);
            }
            if (sm == 0.0) {
                jacobiSort(d, &v);
                return;
            }
            if (i < 4)
                tresh = 0.2 * sm / (n * n);
            else
                tresh = 0.0;

            //calcula la suma de la magnitud de la diagonal
            for (ip = 0; ip < n - 1; ip++) {
                for (iq = ip + 1; iq < n; iq++) {
                    g = 100.0 * abs(a[ip][iq]);
                    if (i > 4 && g <= EPS * abs(d[ip]) && g <= EPS * abs(d[iq]))
                        a[ip][iq] = 0.0;
                    else if (abs(a[ip][iq]) > tresh) {
                        h = d[iq] - d[ip];
                        if (g <= EPS * abs(h))
                            t = (a[ip][iq]) / h;
                        else {
                            theta = 0.5 * h / (a[ip][iq]);
                            t = 1.0 / (abs(theta) + sqrt(1.0 + theta * theta));
                            if (theta < 0.0) t = -t;
                        }
                        c = 1.0 / sqrt(1 + t * t);
                        s = t * c;
                        tau = s / (1.0 + c);
                        h = t * a[ip][iq];
                        z[ip] -= h;
                        z[iq] += h;
                        d[ip] -= h;
                        d[iq] += h;
                        a[ip][iq] = 0.0;

                        // decide si realizar las rotaciones
                        for (unsigned j = 0; j < ip; j++)
                            rot(a, s, tau, j, ip, j, iq);
                        for (unsigned j = ip + 1; j < iq; j++)
                            rot(a, s, tau, ip, j, j, iq);
                        for (unsigned j = iq + 1; j < n; j++)
                            rot(a, s, tau, ip, j, iq, j);
                        for (unsigned j = 0; j < n; j++)
                            rot(v, s, tau, j, ip, j, iq);
                        ++nrot;
                    }
                }
            }
            for (ip = 0; ip < n; ip++) {
                b[ip] += z[ip];
                d[ip] = b[ip];
                z[ip] = 0.0;
            }

        }
        throw ("Limite de interaciones alcanzado en el metodo de Jacobi");
    }

    /**
     * Se encarga de realizar las rotaciones
     * @param a
     * @param s
     * @param tau
     * @param i
     * @param j
     * @param k
     * @param l
     */
    inline void rot(anpi::Matrix<double> &a, const double s, const double tau, const int i,
                    const int j, const int k, const int l) {
        double g = a[i][j];
        double h = a[k][l];
        a[i][j] = g - s * (h + g * tau);
        a[k][l] = h + s * (g - h * tau);
    }
};

/**
 * Funcion encargada de llamar a Jacobi para la resolucion de eigensistemas. Retorna los eigenvectores y eigenvalores
 * @tparam T
 * @param A
 * @param val
 * @param E
 */
template<typename T>
void jacobi(anpi::Matrix <T> &A, std::vector <T> &val, anpi::Matrix <T> &E) {

    JacobiStruct jac(A);
    E = jac.v;
    val = jac.d;

}

/**
 * Funcion encargada calcular la norma de la resta de dos matrices
 * @tparam T
 * @param A
 * @param B
 * @return
 */
template<typename T>
T normaMatrix(anpi::Matrix <T> A, anpi::Matrix <T> B) {
    anpi::Matrix <T> temp = A - B;
    int size = A.rows();
    T norm = T(0);

    for (unsigned int i = 0; i < size; i++) {
        for (unsigned int j = 0; j < size; j++) {
            norm += temp[i][j] * temp[i][j];
        }
    }

    norm = sqrt(norm);
    return norm;
}

/**
 * Funcion que retorna la matriz original reconstruida con los eigenvectores y eigenvalores
 * @tparam T
 * @param E
 * @param D
 * @return
 */
template<typename T>
anpi::Matrix <T> reconstruirMatrix(anpi::Matrix <T> E, std::vector <T> D) {
    unsigned int size = E.rows();
    anpi::Matrix<T> A(size, size, anpi::Matrix<T>::DoNotInitialize);
    anpi::Matrix<T> MD(size, size, anpi::Matrix<T>::DoNotInitialize);
    anpi::Matrix<T> ET(size, size, anpi::Matrix<T>::DoNotInitialize);

    for (unsigned int k = 0; k < size; k++) {
        MD[k][k] = D[k];
    }

    for (unsigned int i = 0; i < size; i++) {
        for (unsigned int j = 0; j < size; j++) {
            ET[j][i] = E[i][j];
        }
    }

    A = E * MD * ET;
    return A;
}

/**
 * Ordena los eigenvalores de forma descendiente y ordena los eigenvectores
 * @tparam T
 * @param val
 * @param E
 */
template<typename T>
void sort(std::vector<T>& val, anpi::Matrix<T>& E){
    unsigned int size = E.rows();
    std::vector<T> valtemp;
    anpi::Matrix<T> Etemp(size, size,anpi:: Matrix<T>::DoNotInitialize);

    for(unsigned int i = 0; i<size; i++){
        valtemp.push_back(val[size-1-i]);
    }

    for(unsigned int i = 0; i<size; i++){
        for(int j = 0; j<size; j++){
            Etemp[i][j] = E[size-1-i][j];
        }
    }

    val = valtemp;
    E = Etemp;

}