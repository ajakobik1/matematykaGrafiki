#include <iostream>
#include <cmath>
#include <vector>
#include <limits>

using namespace std;

// Structure for 3D vector
struct Vec3 {
    double x, y, z;

    Vec3(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}

    Vec3 operator+(const Vec3& v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
    Vec3 operator-(const Vec3& v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
    Vec3 operator*(double t) const { return Vec3(x * t, y * t, z * t); }
    Vec3 operator/(double t) const { return Vec3(x / t, y / t, z / t); }

    double length() const { return sqrt(x*x + y*y + z*z); }
    Vec3 normalize() const {
        double len = length();
        return len > 0 ? (*this) / len : Vec3(0, 0, 0);
    }

    double dot(const Vec3& v) const { return x*v.x + y*v.y + z*v.z; }
};

// Struktura reprezentująca promień (jak promień światła lub laser)
struct Ray {
    Vec3 origin;    // Punkt początkowy (miejsce, z którego wystrzeliwujemy promień, np. oko kamery)
    Vec3 direction; // Kierunek, w którym leci promień (wektor)

    // Konstruktor promienia:
    // Pobiera punkt startowy (o) i kierunek (d).
    // Kierunek jest od razu "normalizowany" (skracany/wydłużany do długości 1),
    // aby ułatwić późniejsze obliczenia odległości.
    Ray(const Vec3& o, const Vec3& d) : origin(o), direction(d.normalize()) {}

    // Funkcja obliczająca punkt na linii promienia w danej odległości 't'
    // t to parametr (liczba), który mówi nam, jak daleko poleciał promień.
    // Wzór: Punkt = Start + (Kierunek * Odległość)
    Vec3 pointAt(double t) const {
        return origin + direction * t;
    }
};

// Cube class (AABB - Axis Aligned Bounding Box)
class Cube {
private:
    Vec3 minPoint; // Punkt o minimalnych współrzędnych (x_min, y_min, z_min)
    Vec3 maxPoint; // Punkt o maksymalnych współrzędnych (x_max, y_max, z_max)

public:
    /**
     * Konstruktor inicjalizujący granice bryły na podstawie punktu centralnego i długości boku.
     */
    Cube(const Vec3& center, double size) {
        double half = size / 2.0;
        minPoint = Vec3(center.x - half, center.y - half, center.z - half);
        maxPoint = Vec3(center.x + half, center.y + half, center.z + half);
    }

    /**
     * Sprawdza przecięcie promienia z trzema parami płaszczyzn ograniczających bryłę (slabs).
     * * @param ray Referencja do obiektu promienia.
     * @param tMin Parametr wyjściowy zwracający najbliższą odległość do punktu trafienia.
     * @return true, jeśli promień przecina objętość bryły; w przeciwnym razie false.
     */
    bool intersect(const Ray& ray, double& tMin) const {
        // Wyznaczenie przedziałów wejścia/wyjścia dla każdej z osi X, Y, Z
        // Wykorzystujemy parametryczne równanie prostej: P(t) = O + tD

        // Analiza dla osi X
        double t1 = (minPoint.x - ray.origin.x) / ray.direction.x;
        double t2 = (maxPoint.x - ray.origin.x) / ray.direction.x;
        double tmin = min(t1, t2);
        double tmax = max(t1, t2);

        // Analiza dla osi Y - wyznaczenie części wspólnej (iloczynu zbiorów) przedziałów t
        t1 = (minPoint.y - ray.origin.y) / ray.direction.y;
        t2 = (maxPoint.y - ray.origin.y) / ray.direction.y;
        tmin = max(tmin, min(t1, t2));
        tmax = min(tmax, max(t1, t2));

        // Analiza dla osi Z
        t1 = (minPoint.z - ray.origin.z) / ray.direction.z;
        t2 = (maxPoint.z - ray.origin.z) / ray.direction.z;
        tmin = max(tmin, min(t1, t2));
        tmax = min(tmax, max(t1, t2));

        // Warunek trafienia: przedział [tmin, tmax] musi być niepusty (tmax >= tmin)
        // oraz znajdować się przed kamerą (tmax >= 0).
        if (tmax >= tmin && tmax >= 0) {
            // Wybór najbliższego punktu trafienia (z uwzględnieniem sytuacji, gdy kamera jest wewnątrz bryły)
            tMin = (tmin >= 0) ? tmin : tmax;
            return true;
        }
        return false;
    }
};

// Camera class
class Camera {
private:
    Vec3 position;
    Vec3 target;
    Vec3 up;
    double fov;

    // Orbit parameters
    double radius;
    double angleH;  // horizontal angle (azimuth)
    double angleV;  // vertical angle (elevation)

public:
    Camera(double r = 5.0, double ah = 0.0, double av = 0.0)
            : radius(r), angleH(ah), angleV(av), target(0, 0, 0), up(0, 1, 0), fov(60.0) {
        updatePosition();
    }

    void updatePosition() {
        // Convert spherical to Cartesian coordinates
        position.x = radius * cos(angleV) * cos(angleH);
        position.y = radius * sin(angleV);
        position.z = radius * cos(angleV) * sin(angleH);
    }

    void rotateHorizontal(double angle) {
        angleH += angle;
        updatePosition();
    }

    void rotateVertical(double angle) {
        angleV += angle;
        updatePosition();
    }

    void setAngleHorizontal(double degrees) {
        angleH = degrees * M_PI / 180.0;
        updatePosition();
    }

    void setAngleVertical(double degrees) {
        angleV = degrees * M_PI / 180.0;
        updatePosition();
    }

    void setZoom(double dist) {
        radius = dist;
        if (radius < 1.0) radius = 1.0;
        if (radius > 20.0) radius = 20.0;
        updatePosition();
    }

    void zoom(double delta) {
        radius += delta;
        if (radius < 1.0) radius = 1.0;
        if (radius > 20.0) radius = 20.0;
        updatePosition();
    }

    /**
   * Generuje promień dla danego punktu na płaszczyźnie obrazu.
   * @param u Znormalizowana współrzędna pozioma [0, 1].
   * @param v Znormalizowana współrzędna pionowa [0, 1].
   */
    Ray getRay(double u, double v, int width, int height) const {
        // 1. Wyznaczenie lokalnego układu współrzędnych kamery (Camera Basis)
        Vec3 forward = (target - position).normalize(); // Oś Z kamery
        Vec3 right = Vec3(forward.z, 0, -forward.x).normalize(); // Oś X kamery (prostopadła do forward i pionu)

        // Iloczyn wektorowy wyznaczający oś Y kamery (górę)
        Vec3 newUp = Vec3(
                forward.y * right.z - forward.z * right.y,
                forward.z * right.x - forward.x * right.z,
                forward.x * right.y - forward.y * right.x
        ).normalize();

        // 2. Obliczenie wymiarów rzutni (Viewport) na podstawie kąta widzenia (FOV)
        double aspectRatio = (double)width / height;
        // Zamiana FOV na radiany i wyznaczenie wysokości płaszczyzny rzutowania
        double viewportHeight = 2.0 * tan(fov * M_PI / 360.0);
        double viewportWidth = viewportHeight * aspectRatio;

        // 3. Mapowanie współrzędnych (u, v) na fizyczne wymiary rzutni
        // Przesunięcie o -0,5 centruje układ współrzędnych na środku ekranu
        double x = (u - 0.5) * viewportWidth;
        double y = (v - 0.5) * viewportHeight;

        // 4. Konstrukcja kierunku promienia jako kombinacji liniowej wektorów bazy
        Vec3 rayDir = forward + right * x + newUp * y;
        return Ray(position, rayDir);
    }

    Vec3 getPosition() const { return position; }
    double getRadius() const { return radius; }
    double getAngleH() const { return angleH * 180.0 / M_PI; }
    double getAngleV() const { return angleV * 180.0 / M_PI; }
};

/**
 * Główna funkcja potoku renderującego (Rendering Pipeline).
 * Przetwarza scenę metodą Ray Castingu, mapując przestrzeń 3D na znaki ASCII.
 */
void render(const Camera& camera, const Cube& cube, int width, int height) {
    // Inicjalizacja bufora ekranu domyślnym znakiem tła
    vector<vector<char>> screen(height, vector<char>(width, '.'));

    // Iteracja po wszystkich pikselach płaszczyzny obrazu
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            // Mapowanie współrzędnych rastrowych (pikseli) na współrzędne znormalizowane [0, 1]
            double u = (double)x / (width - 1);
            double v = (double)(height - 1 - y) / (height - 1); // Korekta orientacji osi Y

            // Generowanie promienia rzutującego dla danego piksela
            Ray ray = camera.getRay(u, v, width, height);

            // Test przecięcia promienia z geometrią sceny (kostką)
            double t;
            if (cube.intersect(ray, t)) {
                // Jeśli nastąpiła kolizja, wpisujemy znak obiektu do bufora
                screen[y][x] = '0';
            }
        }
    }

    // Wyświetlanie zawartości bufora w konsoli z ramką pomocniczą
    system("clear || cls"); // Odświeżanie ekranu

    // Rysowanie górnej krawędzi ramki
    cout << "+";
    for (int i = 0; i < width; i++) cout << "-";
    cout << "+\n";

    // Wypisywanie wierszy bufora ekranu
    for (int y = 0; y < height; y++) {
        cout << "|";
        for (int x = 0; x < width; x++) {
            cout << screen[y][x];
        }
        cout << "|\n";
    }

    // Rysowanie dolnej krawędzi oraz interfejsu użytkownika
    cout << "+";
    for (int i = 0; i < width; i++) cout << "-";
    cout << "+\n";

    // Wyświetlanie metadanych kamery (parametry sferyczne)
    cout << "\nCamera: distance=" << camera.getRadius()
         << " | horizontal=" << (int)round(camera.getAngleH()) << " deg"
         << " | vertical=" << (int)round(camera.getAngleV()) << " deg\n";
    cout << "\nCommands:\n";
    cout << "  Single key:  w/s (vertical) | a/d (horizontal) | +/- (zoom)\n";
    cout << "  Set angle:   h <degrees> (horizontal) | v <degrees> (vertical)\n";
    cout << "  Set zoom:    z <distance>\n";
    cout << "  Combined:    r <h_deg> <v_deg> (rotate both)\n";
    cout << "  Quit:        q\n";
    cout << "\nCommand: ";

}

int main() {
    const int WIDTH = 60;
    const int HEIGHT = 60;

    // Creating cube at the center of coordinate system
    Cube cube(Vec3(0, 0, 0), 2.0);

    // Creating camera - positioned to see cube from the side
    Camera camera(5.0, M_PI / 2.0, 0.0);

    string cmd;
    bool running = true;

    while (running) {
        render(camera, cube, WIDTH, HEIGHT);

        cin >> cmd;

        if (cmd == "w") {
            camera.rotateVertical(5.0 * M_PI / 180.0);
        } else if (cmd == "s") {
            camera.rotateVertical(-5.0 * M_PI / 180.0);
        } else if (cmd == "a") {
            camera.rotateHorizontal(-5.0 * M_PI / 180.0);
        } else if (cmd == "d") {
            camera.rotateHorizontal(5.0 * M_PI / 180.0);
        } else if (cmd == "+") {
            camera.zoom(-0.5);
        } else if (cmd == "-") {
            camera.zoom(0.5);
        } else if (cmd == "h") {
            double angle;
            cin >> angle;
            camera.setAngleHorizontal(angle);
        } else if (cmd == "v") {
            double angle;
            cin >> angle;
            camera.setAngleVertical(angle);
        } else if (cmd == "z") {
            double dist;
            cin >> dist;
            camera.setZoom(dist);
        } else if (cmd == "r") {
            double h_angle, v_angle;
            cin >> h_angle >> v_angle;
            camera.setAngleHorizontal(h_angle);
            camera.setAngleVertical(v_angle);
        } else if (cmd == "q") {
            running = false;
        } else {
            cout << "Unknown command!\n";
            cin.ignore(10000, '\n'); // Clear input buffer
        }
    }

    cout << "\nGoodbye!\n";
    return 0;
}