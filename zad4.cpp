#include <iostream>
#include <fstream> // Biblioteka do obslugi plikow
#include <cmath>
#include <vector>
#include <string>
#include <iomanip>

using namespace std;

// Stała dla liczby PI
const double PI = 3.14159265358979323846;

// --- KLASA WEKTORA I FUNKCJE POMOCNICZE ---

struct Vector3 {
    double x, y, z;

    Vector3(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}

    // Operacje podstawowe
    Vector3 operator+(const Vector3& other) const { return Vector3(x + other.x, y + other.y, z + other.z); }
    Vector3 operator-(const Vector3& other) const { return Vector3(x - other.x, y - other.y, z - other.z); }
    Vector3 operator*(double scalar) const { return Vector3(x * scalar, y * scalar, z * scalar); }
    Vector3 operator-() const { return Vector3(-x, -y, -z); }

    // Długość wektora
    double length() const { return sqrt(x * x + y * y + z * z); }

    // Normalizacja
    Vector3 normalize() const {
        double len = length();
        if (len == 0) return *this;
        return Vector3(x / len, y / len, z / len);
    }

    // Wypisywanie do dowolnego strumienia (plik lub konsola)
    friend ostream& operator<<(ostream& os, const Vector3& v) {
        os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
        return os;
    }
};

// Iloczyn skalarny (Dot Product)
double dotProduct(const Vector3& a, const Vector3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// Iloczyn wektorowy (Cross Product)
Vector3 crossProduct(const Vector3& a, const Vector3& b) {
    return Vector3(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}

// Konwersja radianów na stopnie
double toDegrees(double radians) {
    return radians * 180.0 / PI;
}

// Struktura Prostej (punkt + wektor kierunkowy)
struct Line {
    Vector3 p; // Punkt początkowy
    Vector3 v; // Wektor kierunkowy

    Vector3 getPointAt(double t) const {
        return p + v * t;
    }
};

// Struktura Płaszczyzny (równanie ogólne Ax + By + Cz + D = 0)
// Normalna n = [A, B, C]
struct Plane {
    Vector3 n; // Normalna
    double d;  // Współczynnik D

    Plane(double a, double b, double c, double d_val) : n(a, b, c), d(d_val) {}
};

// --- FUNKCJE ROZWIĄZUJĄCE ZADANIA ---

// Zadanie 1 & 7: Przecięcie prostych (lub odcinków)
bool intersectLines(const Line& l1, const Line& l2, Vector3& intersectionPoint, bool checkSegments = false) {
    Vector3 p1p2 = l2.p - l1.p;
    Vector3 v1xv2 = crossProduct(l1.v, l2.v);

    if (v1xv2.length() < 1e-9) return false; // Równoległe

    double det = v1xv2.length() * v1xv2.length();
    double t = dotProduct(crossProduct(p1p2, l2.v), v1xv2) / det;
    double u = dotProduct(crossProduct(p1p2, l1.v), v1xv2) / det;

    if (checkSegments) {
        if (t < 0 || t > 1 || u < 0 || u > 1) return false;
    }

    intersectionPoint = l1.getPointAt(t);
    return true;
}

// Zadanie 2: Kąt między prostymi
double angleBetweenLines(const Line& l1, const Line& l2) {
    double dot = dotProduct(l1.v, l2.v);
    double lenMul = l1.v.length() * l2.v.length();
    return toDegrees(acos(abs(dot) / lenMul));
}

// Zadanie 3: Przecięcie prostej i płaszczyzny
bool intersectLinePlane(const Line& l, const Plane& pl, Vector3& intersectionPoint) {
    double denominator = dotProduct(pl.n, l.v);
    if (abs(denominator) < 1e-9) return false;

    double t = -(dotProduct(pl.n, l.p) + pl.d) / denominator;
    intersectionPoint = l.getPointAt(t);
    return true;
}

// Zadanie 4: Kąt między prostą a płaszczyzną
double angleLinePlane(const Line& l, const Plane& pl) {
    double nominator = abs(dotProduct(pl.n, l.v));
    double denominator = pl.n.length() * l.v.length();
    return toDegrees(asin(nominator / denominator));
}

// Zadanie 5: Prosta przecięcia płaszczyzn
bool intersectPlanes(const Plane& p1, const Plane& p2, Line& resultLine) {
    Vector3 dir = crossProduct(p1.n, p2.n);
    if (dir.length() < 1e-9) return false;

    resultLine.v = dir;

    double det = p1.n.x * p2.n.y - p1.n.y * p2.n.x;

    if (abs(det) > 1e-9) {
        resultLine.p.z = 0;
        resultLine.p.x = (p1.n.y * p2.d - p2.n.y * p1.d) / det;
        resultLine.p.y = (p2.n.x * p1.d - p1.n.x * p2.d) / det;
    } else {
        det = p1.n.x * p2.n.z - p1.n.z * p2.n.x;
        if(abs(det) > 1e-9) {
             resultLine.p.y = 0;
             resultLine.p.x = (p1.n.z * p2.d - p2.n.z * p1.d) / det;
             resultLine.p.z = (p2.n.x * p1.d - p1.n.x * p2.d) / det;
        } else {
            det = p1.n.y * p2.n.z - p1.n.z * p2.n.y;
             resultLine.p.x = 0;
             resultLine.p.y = (p1.n.z * p2.d - p2.n.z * p1.d) / det;
             resultLine.p.z = (p2.n.y * p1.d - p1.n.y * p2.d) / det;
        }
    }
    return true;
}

// Zadanie 6: Kąt między płaszczyznami
double angleBetweenPlanes(const Plane& p1, const Plane& p2) {
    double dot = dotProduct(p1.n, p2.n);
    double lenMul = p1.n.length() * p2.n.length();
    return toDegrees(acos(abs(dot) / lenMul));
}

// Strona 2, Zadanie 1: Przecięcie sfery i prostej
// Funkcja przyjmuje teraz strumień wyjściowy (plik), żeby zapisać do niego wynik
void intersectSphereLine(const Vector3& center, double r, const Line& l, ostream& out) {
    Vector3 L = l.p - center;
    double a = dotProduct(l.v, l.v);
    double b = 2 * dotProduct(L, l.v);
    double c = dotProduct(L, L) - (r * r);

    double delta = b * b - 4 * a * c;

    out << "ZADANIE (str. 2) nr 1: Punkt przeciecia sfery i prostej" << endl;
    out << "DANE: Sfera r=" << r << ", Prosta przez A i A'" << endl;

    if (delta < 0) {
        out << "WYNIK: Brak punktow przeciecia (delta < 0)." << endl;
    } else if (abs(delta) < 1e-9) {
        double t = -b / (2 * a);
        out << "WYNIK: Jeden punkt (styczny): " << l.getPointAt(t) << endl;
    } else {
        double t1 = (-b - sqrt(delta)) / (2 * a);
        double t2 = (-b + sqrt(delta)) / (2 * a);
        out << "WYNIK: Dwa punkty przeciecia:" << endl;
        out << "       Punkt 1: " << l.getPointAt(t1) << endl;
        out << "       Punkt 2: " << l.getPointAt(t2) << endl;
    }
    out << "-----------------------------------" << endl;
}

// --- MAIN ---

int main() {
    // Otwarcie pliku do zapisu
    ofstream resultFile("wyniki.txt");

    if (!resultFile.is_open()) {
        cerr << "Blad: Nie udalo sie otworzyc pliku wyniki.txt do zapisu!" << endl;
        return 1;
    }

    // Ustawienie formatowania liczb w pliku
    resultFile << fixed << setprecision(2);

    cout << "Obliczanie rozpoczete..." << endl;

    // --- ZADANIE 1: PUNKT PRZECIĘCIA PROSTYCH ---
    Line line1A = { Vector3(-2, 4, 0), Vector3(3, 1, 5) };
    Line line1B = { Vector3(-2, 4, 0), Vector3(1, -5, 3) };

    Vector3 intersection1;
    resultFile << "ZADANIE 1: Znajdz punkt przeciecia prostych" << endl;
    resultFile << "Prosta A: (x+2)/3 = y-4 = z/5" << endl;
    resultFile << "Prosta B: (x+2)/1 = -1/5(y-4) = z/3" << endl;
    if (intersectLines(line1A, line1B, intersection1)) {
        resultFile << "WYNIK: Punkt przeciecia to " << intersection1 << endl;
    } else {
        resultFile << "WYNIK: Proste sie nie przecinaja." << endl;
    }
    resultFile << "-----------------------------------" << endl;

    // --- ZADANIE 2: KĄT MIĘDZY PROSTYMI ---
    resultFile << "ZADANIE 2: Znajdz kat miedzy prostymi z zadania 1" << endl;
    resultFile << "WYNIK: " << angleBetweenLines(line1A, line1B) << " stopni" << endl;
    resultFile << "-----------------------------------" << endl;

    // --- ZADANIE 3: PRZECIĘCIE PROSTEJ I PŁASZCZYZNY ---
    Line line3 = { Vector3(-2, 2, -1), Vector3(3, -1, 2) };
    Plane plane3(2, 3, 3, -8);

    Vector3 intersection3;
    resultFile << "ZADANIE 3: Znajdz punkt przeciecia prostej i plaszczyzny" << endl;
    resultFile << "Prosta: parametryczna t, Płaszczyzna: 2x+3y+3z-8=0" << endl;
    if (intersectLinePlane(line3, plane3, intersection3)) {
        resultFile << "WYNIK: Punkt przeciecia to " << intersection3 << endl;
    } else {
        resultFile << "WYNIK: Brak (rownolegle)." << endl;
    }
    resultFile << "-----------------------------------" << endl;

    // --- ZADANIE 4: KĄT MIĘDZY PROSTĄ A PŁASZCZYZNĄ ---
    resultFile << "ZADANIE 4: Znajdz kat pomiedzy prosta a plaszczyzna z zadania 2 (wg treści zadania 3)" << endl;
    resultFile << "WYNIK: " << angleLinePlane(line3, plane3) << " stopni" << endl;
    resultFile << "-----------------------------------" << endl;

    // --- ZADANIE 5: PROSTA PRZECIĘCIA PŁASZCZYZN ---
    Plane plane5A(2, -1, 1, -8);
    Plane plane5B(4, 3, 1, 14);

    Line intersectionLine5;
    resultFile << "ZADANIE 5: Znajdz prosta przeciecia plaszczyzn" << endl;
    resultFile << "Plaszczyzna A: 2x-y+z-8=0" << endl;
    resultFile << "Plaszczyzna B: 4x+3y+z+14=0" << endl;
    if (intersectPlanes(plane5A, plane5B, intersectionLine5)) {
        resultFile << "WYNIK: Rownanie parametryczne prostej:" << endl;
        resultFile << "       Punkt P" << intersectionLine5.p << " + t * Wektor V" << intersectionLine5.v << endl;
    } else {
        resultFile << "WYNIK: Plaszczyzny sa rownolegle." << endl;
    }
    resultFile << "-----------------------------------" << endl;

    // --- ZADANIE 6: KĄT MIĘDZY PŁASZCZYZNAMI ---
    resultFile << "ZADANIE 6: Znajdz kat pomiedzy plaszczyznami z zadania 5" << endl;
    resultFile << "WYNIK: " << angleBetweenPlanes(plane5A, plane5B) << " stopni" << endl;
    resultFile << "-----------------------------------" << endl;

    // --- ZADANIE 7: PRZECIĘCIE ODCINKÓW ---
    Vector3 startA(5,5,4), endA(10,10,6);
    Vector3 startB(5,5,5), endB(10,10,3);
    Line segA = { startA, endA - startA };
    Line segB = { startB, endB - startB };

    Vector3 intersection7;
    resultFile << "ZADANIE 7: Znajdz punkt przeciecia dwoch odcinkow" << endl;
    resultFile << "Odcinek A: (5,5,4)->(10,10,6)" << endl;
    resultFile << "Odcinek B: (5,5,5)->(10,10,3)" << endl;
    if (intersectLines(segA, segB, intersection7, true)) {
        resultFile << "WYNIK: Punkt przeciecia to " << intersection7 << endl;
    } else {
        resultFile << "WYNIK: Odcinki sie nie przecinaja (punkt wspolny poza zakresem odcinka)." << endl;
    }
    resultFile << "-----------------------------------" << endl;

    // --- STRONA 2, ZADANIE 1: SFERA I PROSTA ---
    Vector3 sphereCenter(0, 0, 0);
    double radius = sqrt(26);
    Vector3 pA(3, -1, -2);
    Vector3 pApr(5, 3, -4);
    Line lineSphere = { pA, pApr - pA };

    // Tutaj przekazujemy obiekt resultFile zamiast wypisywac w funkcji na cout
    intersectSphereLine(sphereCenter, radius, lineSphere, resultFile);

    // Zamkniecie pliku
    resultFile.close();

    cout << "Gotowe! Wyniki zapisano w pliku wyniki.txt" << endl;

    return 0;
}