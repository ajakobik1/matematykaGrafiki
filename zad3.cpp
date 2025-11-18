#include <iostream>
#include <fstream>
#include <cmath>

class Quaternion {
public:
    double w, x, y, z;

    Quaternion(double w=0, double x=0, double y=0, double z=0)
        : w(w), x(x), y(y), z(z) {}

    // Dodawanie
    Quaternion operator+(const Quaternion &q) const {
        return Quaternion(w + q.w, x + q.x, y + q.y, z + q.z);
    }

    // Odejmowanie
    Quaternion operator-(const Quaternion &q) const {
        return Quaternion(w - q.w, x - q.x, y - q.y, z - q.z);
    }

    // Mnożenie
    Quaternion operator*(const Quaternion &q) const {
        return Quaternion(
            w*q.w - x*q.x - y*q.y - z*q.z,
            w*q.x + x*q.w + y*q.z - z*q.y,
            w*q.y - x*q.z + y*q.w + z*q.x,
            w*q.z + x*q.y - y*q.x + z*q.w
        );
    }

    // Sprzężenie
    Quaternion conjugate() const {
        return Quaternion(w, -x, -y, -z);
    }

    // Norma
    double norm() const {
        return std::sqrt(w*w + x*x + y*y + z*z);
    }

    // Normalizacja
    Quaternion normalized() const {
        double n = norm();
        return Quaternion(w/n, x/n, y/n, z/n);
    }

    // Zapis do streamu
    void write(std::ostream &os, std::string name = "") const {
        os << name << "(" << w << ", " << x << ", " << y << ", " << z << ")\n";
    }
};

// Obrót punktu (wektora) przez kwaternion
Quaternion rotatePoint(const Quaternion &q, const Quaternion &p) {
    Quaternion qn = q.normalized();
    return qn * p * qn.conjugate();
}

int main() {
    std::ofstream file("wyniki.txt");
    if (!file.is_open()) {
        std::cerr << "Nie udalo sie otworzyc pliku!\n";
        return 1;
    }

    file << "=== TESTY KWATERNIONOW ===\n\n";

    Quaternion a(1, 2, 3, 4);
    Quaternion b(0.5, -1, 2, 0);

    file << "--- A i B ---\n";
    a.write(file, "A = ");
    b.write(file, "B = ");
    file << "\n";

    // Testy działań
    file << "--- Dzialania ---\n";
    (a + b).write(file, "A + B = ");
    (a - b).write(file, "A - B = ");
    (a * b).write(file, "A * B = ");
    (b * a).write(file, "B * A = ");

    file << "\n--- Sprzezenie ---\n";
    a.conjugate().write(file, "conj(A) = ");
    b.conjugate().write(file, "conj(B) = ");

    file << "\n--- Norma ---\n";
    file << "||A|| = " << a.norm() << "\n";
    file << "||B|| = " << b.norm() << "\n";

    file << "\n--- Normalizacja ---\n";
    a.normalized().write(file, "A_norm = ");
    b.normalized().write(file, "B_norm = ");

    // Niezmienność mnożenia
    file << "\n=== BRAK PRZEMIENNOSCI ===\n";
    Quaternion ab = a * b;
    Quaternion ba = b * a;

    ab.write(file, "A * B = ");
    ba.write(file, "B * A = ");
    file << "Czy A * B = B * A?  ->  NIE\n";

    // Obrót punktu
    file << "\n=== OBRÓT PUNKTU ===\n";
    file << "Punkt: (-1, -1, -1)\n";
    file << "Obrot: 270 stopni wokol osi X\n";

    double angle = 270 * M_PI / 180.0;
    double half = angle / 2.0;

    // Oś X = (1,0,0)
    Quaternion q_rot(std::cos(half), std::sin(half), 0, 0);
    Quaternion point(0, -1, -1, -1);

    Quaternion rotated = rotatePoint(q_rot, point);

    rotated.write(file, "Po obrocie punkt = ");

    file.close();
    std::cout << "Wyniki zapisane w pliku wyniki.txt\n";

    return 0;
}
