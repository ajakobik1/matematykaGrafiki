#include <iostream>
#include <cmath>

class Vector3 {
public:
    float x, y, z;

    // konstruktory
    Vector3() : x(0), y(0), z(0) {}
    Vector3(float x, float y, float z) : x(x), y(y), z(z) {}

    // dodawanie i odejmowanie
    Vector3 operator+(const Vector3& v) const {
        return Vector3(x + v.x, y + v.y, z + v.z);
    }

    Vector3 operator-(const Vector3& v) const {
        return Vector3(x - v.x, y - v.y, z - v.z);
    }

    // mnożenie i dzielenie przez skalar
    Vector3 operator*(float scalar) const {
        return Vector3(x * scalar, y * scalar, z * scalar);
    }

    Vector3 operator/(float scalar) const {
        return Vector3(x / scalar, y / scalar, z / scalar);
    }

    // iloczyn skalarny
    float dot(const Vector3& v) const {
        return x * v.x + y * v.y + z * v.z;
    }

    // iloczyn wektorowy
    Vector3 cross(const Vector3& v) const {
        return Vector3(
            y * v.z - z * v.y,
            z * v.x - x * v.z,
            x * v.y - y * v.x
        );
    }

    // długość wektora
    float length() const {
        return std::sqrt(x * x + y * y + z * z);
    }

    // normalizacja wektora
    Vector3 normalize() const {
        float len = length();
        if (len == 0) return Vector3(0, 0, 0);
        return *this / len;
    }

    // wyświetlanie wektora
    void print(const std::string& name = "") const {
        if (!name.empty())
            std::cout << name << " = ";
        std::cout << "[" << x << ", " << y << ", " << z << "]" << std::endl;
    }
};

float calculateAngle(const Vector3 &v1, const Vector3 &v2) {
    float dotProd = v1.dot(v2);
    float len1 = v1.length();
    float len2 = v2.length();
    float cosTheta = dotProd / (len1 * len2);
    float angleRad = std::acos(cosTheta);
    float angleDeg = angleRad * 180.0f / M_PI;
    return angleDeg;
}

int main() {
    // 1. czy dodawanie jest przemienne?
    Vector3 a(1, 2, 3);
    Vector3 b(4, 5, 6);
    Vector3 sum1 = a + b;
    Vector3 sum2 = b + a;

    sum1.print("a + b");
    sum2.print("b + a");

    std::cout << "Czy a + b = b + a ? "
              << ((sum1.x == sum2.x && sum1.y == sum2.y && sum1.z == sum2.z) ? "TAK" : "NIE")
              << std::endl;

    std::cout << "-------------------------------------------------\n";

    // 2. Kat pomiedzy [0,3,0] a [5,5,0]
    Vector3 v1(0, 3, 0);
    Vector3 v2(5, 5, 0);

    float angleDeg = calculateAngle(v1, v2);

    std::cout << "Kat pomiedzy [0,3,0] a [5,5,0] = " << angleDeg << " stopni" << std::endl;

    std::cout << "-------------------------------------------------\n";

    // 3. Wektor prostopadly do [4,5,1] i [4,1,3]
    Vector3 v3(4, 5, 1);
    Vector3 v4(4, 1, 3);
    Vector3 perpendicular = v3.cross(v4);

    perpendicular.print("Wektor prostopadly");

    std::cout << "-------------------------------------------------\n";

    // 4. Normalizacja powstalego wektora
    Vector3 normalized = perpendicular.normalize();
    normalized.print("Wektor znormalizowany");

    std::cout << "Dlugosc wektora znormalizowanego = " << normalized.length() << std::endl;

    return 0;
}


