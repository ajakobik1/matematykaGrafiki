#include <iostream>
#include <cmath>
#include <iomanip>

#define DEG2RAD(x) ((x) * M_PI / 180.0)

class Macierz {
public:
    double m[4][4];

    // Konstruktor domyslny - macierz zerowa
    Macierz() {
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                m[i][j] = 0.0;
    }

    // Macierz jednostkowa
    static Macierz Jednostkowa() {
        Macierz I;
        for (int i = 0; i < 4; ++i)
            I.m[i][i] = 1.0;
        return I;
    }

    // Dodawanie macierzy
    Macierz operator+(const Macierz& B) const {
        Macierz wynik;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                wynik.m[i][j] = m[i][j] + B.m[i][j];
        return wynik;
    }

    // Odejmowanie macierzy
    Macierz operator-(const Macierz& B) const {
        Macierz wynik;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                wynik.m[i][j] = m[i][j] - B.m[i][j];
        return wynik;
    }

    // Mnozenie przez skalar
    Macierz operator*(double s) const {
        Macierz wynik;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                wynik.m[i][j] = m[i][j] * s;
        return wynik;
    }

    // Mnozenie macierzy
    Macierz operator*(const Macierz& B) const {
        Macierz wynik;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                for (int k = 0; k < 4; ++k)
                    wynik.m[i][j] += m[i][k] * B.m[k][j];
        return wynik;
    }

    // Transpozycja
    Macierz Transpozycja() const {
        Macierz wynik;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                wynik.m[i][j] = m[j][i];
        return wynik;
    }

    // Macierz translacji
    static Macierz Translacja(double a, double b, double c) {
        Macierz T = Jednostkowa();
        T.m[0][3] = a;
        T.m[1][3] = b;
        T.m[2][3] = c;
        return T;
    }

    // Macierz skalowania
    static Macierz Skalowanie(double sx, double sy, double sz) {
        Macierz S = Jednostkowa();
        S.m[0][0] = sx;
        S.m[1][1] = sy;
        S.m[2][2] = sz;
        return S;
    }

    // Obrot wokol osi X
    static Macierz RotX(double kat) {
        Macierz R = Jednostkowa();
        double r = DEG2RAD(kat);
        R.m[1][1] = cos(r);
        R.m[1][2] = -sin(r);
        R.m[2][1] = sin(r);
        R.m[2][2] = cos(r);
        return R;
    }

    // Obrot wokol osi Y
    static Macierz RotY(double kat) {
        Macierz R = Jednostkowa();
        double r = DEG2RAD(kat);
        R.m[0][0] = cos(r);
        R.m[0][2] = sin(r);
        R.m[2][0] = -sin(r);
        R.m[2][2] = cos(r);
        return R;
    }

    // Obrot wokol osi Z
    static Macierz RotZ(double kat) {
        Macierz R = Jednostkowa();
        double r = DEG2RAD(kat);
        R.m[0][0] = cos(r);
        R.m[0][1] = -sin(r);
        R.m[1][0] = sin(r);
        R.m[1][1] = cos(r);
        return R;
    }

    // Wyswietlenie macierzy
    void Pokaz() const {
        std::cout << std::fixed << std::setprecision(2);
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j)
                std::cout << std::setw(8) << m[i][j] << " ";
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
};

// Mnozenie wektora przez macierz
void TransformujWektor(double v[4], const Macierz& M) {
    double wynik[4] = {0, 0, 0, 0};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            wynik[i] += M.m[i][j] * v[j];

    std::cout << "Wynik transformacji wektora: [";
    for (int i = 0; i < 4; ++i)
        std::cout << wynik[i] << (i < 3 ? ", " : "]\n");
}

int main() {
    std::cout << "=== Test macierzy ===\n";

    Macierz I = Macierz::Jednostkowa();
    std::cout << "Macierz jednostkowa:\n"; I.Pokaz();

    Macierz A = Macierz::Skalowanie(2, 3, 4);
    Macierz B = Macierz::Translacja(5, 0, -2);

    std::cout << "Macierz A (skalowanie):\n"; A.Pokaz();
    std::cout << "Macierz B (translacja):\n"; B.Pokaz();

    std::cout << "A * B:\n"; (A * B).Pokaz();
    std::cout << "B * A:\n"; (B * A).Pokaz();

    std::cout << "Dowod braku przemiennosci: A*B != B*A\n\n";

    // Obrot wektora [1,0,0,1] o 90 stopni wokol osi Y
    double v[4] = {1, 0, 0, 1};
    std::cout << "Wektor przed obrotem: [1, 0, 0, 1]\n";
    Macierz R = Macierz::RotY(90);
    TransformujWektor(v, R);

    return 0;
}
