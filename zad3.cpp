#include <iostream>
#include <fstream>
#include <cmath>

class Quaternion {
public:
    double w, x, y, z;

    Quaternion(double w=0, double x=0, double y=0, double z=0)
        : w(w), x(x), y(y), z(z) {}

    // PODSTAWOWE DZIAŁANIA

    // Dodawanie
    Quaternion operator+(const Quaternion &q) const {
        return Quaternion(w + q.w, x + q.x, y + q.y, z + q.z);
    }

    // Odejmowanie
    Quaternion operator-(const Quaternion &q) const {
        return Quaternion(w - q.w, x - q.x, y - q.y, z - q.z);
    }

    // Mnożenie (zgodne z notatkami)
    Quaternion operator*(const Quaternion &q) const {
        return Quaternion(
            w*q.w - x*q.x - y*q.y - z*q.z,
            w*q.x + x*q.w + y*q.z - z*q.y,
            w*q.y - x*q.z + y*q.w + z*q.x,
            w*q.z + x*q.y - y*q.x + z*q.w
        );
    }

    // SPRZĘŻENIE (a − v)
    Quaternion conjugate() const {
        return Quaternion(w, -x, -y, -z);
    }

    // NORMA ||q||
    double norm() const {
        return std::sqrt(w*w + x*x + y*y + z*z);
    }

    // Kwadrat normy (potrzebny do odwrotności)
    double normSquared() const {
        return w*w + x*x + y*y + z*z;
    }

    // ODWROTNOŚĆ q⁻¹ = q* / ||q||²
    Quaternion inverse() const {
        double n2 = normSquared();
        Quaternion conj = conjugate();
        return Quaternion(conj.w / n2, conj.x / n2, conj.y / n2, conj.z / n2);
    }

    // DZIELENIE q1 / q2 = q1 * q2⁻¹
    Quaternion operator/(const Quaternion &q) const {
        return (*this) * q.inverse();
    }

    // Zapis kwaternionu do pliku
    void write(std::ostream &os, const std::string &label = "") const {
        os << label << "(" << w << ", " << x << ", " << y << ", " << z << ")\n";
    }
};

// Obrót punktu wektora kwaternionem: v' = q v q*
Quaternion rotatePoint(const Quaternion &q, const Quaternion &p) {
    Quaternion qn = q;
    return qn * p * qn.conjugate();
}

// MAIN — TESTY + ZAPIS DO PLIKU

int main() {
    std::ofstream file("wyniki.txt");
    if (!file.is_open()) {
        std::cerr << "Błąd: nie mogę otworzyć pliku wyniki.txt\n";
        return 1;
    }

    file << "================= TESTY KWATERNIONÓW =================\n\n";

    Quaternion A(1, 2, 3, 4);
    Quaternion B(0.5, -1, 2, 0);

    file << "--- Kwaterniony testowe ---\n";
    A.write(file, "A = ");
    B.write(file, "B = ");
    file << "\n";

    // DZIAŁANIA

    file << "--- Dodawanie ---\n";
    (A + B).write(file, "A + B = ");

    file << "\n--- Odejmowanie ---\n";
    (A - B).write(file, "A - B = ");

    file << "\n--- Mnożenie ---\n";
    (A * B).write(file, "A * B = ");
    (B * A).write(file, "B * A = ");

    file << "\n--- Sprzężenie ---\n";
    A.conjugate().write(file, "conj(A) = ");
    B.conjugate().write(file, "conj(B) = ");

    file << "\n--- Norma i norma^2 ---\n";
    file << "||A|| = " << A.norm() << "\n";
    file << "||B|| = " << B.norm() << "\n";
    file << "||A||^2 = " << A.normSquared() << "\n";
    file << "||B||^2 = " << B.normSquared() << "\n";

    file << "\n--- Odwrotność ---\n";
    A.inverse().write(file, "A^-1 = ");
    B.inverse().write(file, "B^-1 = ");

    file << "\n--- Dzielenie ---\n";
    (A / B).write(file, "A / B = ");
    (B / A).write(file, "B / A = ");

    file << "\nCzy A/B * B = A ? -> Sprawdzamy\n";
    ((A / B) * B).write(file, "(A/B) * B = ");

    // BRAK PRZEMIENNOŚCI
    file << "\n--- Brak przemienności ---\n";
    Quaternion ab = A * B;
    Quaternion ba = B * A;

    ab.write(file, "A * B = ");
    ba.write(file, "B * A = ");
    file << "Czy A * B = B * A?  ->  NIE\n";

    // OBRÓT PUNKTU
    file << "\n--- Obrót punktu (-1, -1, -1) o 270° wokół osi X ---\n";

    double angle = 270 * M_PI / 180.0;
    Quaternion q_rot(std::cos(angle/2), std::sin(angle/2), 0, 0);

    Quaternion point(0, -1, -1, -1);
    Quaternion rotated = rotatePoint(q_rot, point);

    rotated.write(file, "Po obrocie = ");

    file.close();
    std::cout << "Wyniki zapisane w wyniki.txt\n";

    return 0;
}
