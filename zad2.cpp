#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

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

    double Det3x3(double a1, double a2, double a3,
                  double b1, double b2, double b3,
                  double c1, double c2, double c3) const {
        return a1 * (b2 * c3 - b3 * c2)
             - a2 * (b1 * c3 - b3 * c1)
             + a3 * (b1 * c2 - b2 * c1);
    }

    Macierz Odwrotna() const {
        Macierz wynik;
        double det = 0.0;

        // obliczamy wyznacznik 4x4
        det =
            m[0][0] * Det3x3(m[1][1], m[1][2], m[1][3],
                             m[2][1], m[2][2], m[2][3],
                             m[3][1], m[3][2], m[3][3])
          - m[0][1] * Det3x3(m[1][0], m[1][2], m[1][3],
                             m[2][0], m[2][2], m[2][3],
                             m[3][0], m[3][2], m[3][3])
          + m[0][2] * Det3x3(m[1][0], m[1][1], m[1][3],
                             m[2][0], m[2][1], m[2][3],
                             m[3][0], m[3][1], m[3][3])
          - m[0][3] * Det3x3(m[1][0], m[1][1], m[1][2],
                             m[2][0], m[2][1], m[2][2],
                             m[3][0], m[3][1], m[3][2]);

        if (fabs(det) < 1e-9) {
            std::cerr << "Macierz nie ma odwrotnosci (det=0)\n";
            return wynik;
        }

        // Obliczamy macierz dopelnien algebraicznych
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                double sub[3][3];
                int subi = 0;
                for (int ii = 0; ii < 4; ++ii) {
                    if (ii == i) continue;
                    int subj = 0;
                    for (int jj = 0; jj < 4; ++jj) {
                        if (jj == j) continue;
                        sub[subi][subj] = m[ii][jj];
                        subj++;
                    }
                    subi++;
                }
                double subDet = Det3x3(sub[0][0], sub[0][1], sub[0][2],
                                       sub[1][0], sub[1][1], sub[1][2],
                                       sub[2][0], sub[2][1], sub[2][2]);
                wynik.m[j][i] = ((i + j) % 2 == 0 ? 1 : -1) * subDet / det;
            }
        }
        return wynik;
    }

    // Zapis macierzy do strumienia (konsola lub plik)
    void Zapisz(std::ostream& out) const {
        out << std::fixed << std::setprecision(2);
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j)
                out << std::setw(8) << m[i][j] << " ";
            out << "\n";
        }
        out << "\n";
    }
};

// Mnozenie wektora przez macierz
void TransformujWektor(double v[4], const Macierz& M, std::ostream& out) {
    double wynik[4] = {0, 0, 0, 0};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            wynik[i] += M.m[i][j] * v[j];

    out << "Wynik transformacji wektora: [";
    for (int i = 0; i < 4; ++i)
        out << wynik[i] << (i < 3 ? ", " : "]\n");
}

int main() {
    std::ofstream out("wyniki.txt");

    out << "=== Test macierzy ===\n\n";

    Macierz I = Macierz::Jednostkowa();
    out << "Macierz jednostkowa:\n"; I.Zapisz(out);

    Macierz A = Macierz::Skalowanie(2, 3, 4);
    Macierz B = Macierz::Translacja(5, 0, -2);

    out << "Macierz A (skalowanie):\n"; A.Zapisz(out);
    out << "Macierz B (translacja):\n"; B.Zapisz(out);

    out << "A * B:\n"; (A * B).Zapisz(out);
    out << "B * A:\n"; (B * A).Zapisz(out);

    out << "Dowod braku przemiennosci: A*B != B*A\n\n";

    // Odwrotność
    out << "--- Odwrotność macierzy testowej ---\n";
    Macierz C;
    C.m[0][0] = 2; C.m[0][1] = 1; C.m[0][2] = 0; C.m[0][3] = 0;
    C.m[1][0] = 0; C.m[1][1] = 1; C.m[1][2] = 1; C.m[1][3] = 0;
    C.m[2][0] = 1; C.m[2][1] = 0; C.m[2][2] = 1; C.m[2][3] = 0;
    C.m[3][3] = 1;

    out << "Macierz C:\n"; C.Zapisz(out);
    Macierz Cinv = C.Odwrotna();
    out << "Macierz odwrotna C^-1:\n"; Cinv.Zapisz(out);
    out << "Sprawdzenie: C * C^-1 =\n"; (C * Cinv).Zapisz(out);

    // Dodawanie
    out << "--- Dodawanie (C + I) ---\n";
    (C + I).Zapisz(out);

    // Odejmowanie
    out << "--- Odejmowanie (C - I) ---\n";
    (C - I).Zapisz(out);

    // Mnożenie przez skalar
    out << "--- Mnozenie przez skalar (C * 2) ---\n";
    (C * 2).Zapisz(out);

    // Transpozycja
    out << "--- Transpozycja macierzy C ---\n";
    C.Transpozycja().Zapisz(out);

    // Obrot wektora [1,0,0,1] o 90 stopni wokol osi Y
    double v[4] = {1, 0, 0, 1};
    out << "Wektor przed obrotem: [1, 0, 0, 1]\n";
    Macierz R = Macierz::RotY(90);
    TransformujWektor(v, R, out);

    std::cout << "Wyniki zapisane do pliku wyniki.txt\n";


    Macierz G;
    G.m[0][0] = 2; G.m[0][1] = 1; G.m[0][2] = 0; G.m[0][3] = 0;
    G.m[1][0] = 0; G.m[1][1] = 1; G.m[1][2] = 1; G.m[1][3] = 0;
    G.m[2][0] = 1; G.m[2][1] = 0; G.m[2][2] = 1; G.m[2][3] = 0;
    G.m[3][3] = 1;

    G = G.RotX(45);
    G.Zapisz(out);

    G = G.RotX(-45);
    G.Zapisz(out);

    G = G.RotX(-45);
    G.Zapisz(out);

    G = G.Odwrotna();
    G.Zapisz(out);

    out.close();
    return 0;
}
