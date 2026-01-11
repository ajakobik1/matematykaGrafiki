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

// Structure for ray
struct Ray {
    Vec3 origin;
    Vec3 direction;

    Ray(const Vec3& o, const Vec3& d) : origin(o), direction(d.normalize()) {}

    Vec3 pointAt(double t) const { return origin + direction * t; }
};

// Cube class (AABB - Axis Aligned Bounding Box)
class Cube {
private:
    Vec3 minPoint;
    Vec3 maxPoint;

public:
    Cube(const Vec3& center, double size) {
        double half = size / 2.0;
        minPoint = Vec3(center.x - half, center.y - half, center.z - half);
        maxPoint = Vec3(center.x + half, center.y + half, center.z + half);
    }

    // Check ray-cube intersection (slab algorithm)
    bool intersect(const Ray& ray, double& tMin) const {
        double t1 = (minPoint.x - ray.origin.x) / ray.direction.x;
        double t2 = (maxPoint.x - ray.origin.x) / ray.direction.x;

        double tmin = min(t1, t2);
        double tmax = max(t1, t2);

        t1 = (minPoint.y - ray.origin.y) / ray.direction.y;
        t2 = (maxPoint.y - ray.origin.y) / ray.direction.y;

        tmin = max(tmin, min(t1, t2));
        tmax = min(tmax, max(t1, t2));

        t1 = (minPoint.z - ray.origin.z) / ray.direction.z;
        t2 = (maxPoint.z - ray.origin.z) / ray.direction.z;

        tmin = max(tmin, min(t1, t2));
        tmax = min(tmax, max(t1, t2));

        if (tmax >= tmin && tmax >= 0) {
            tMin = tmin >= 0 ? tmin : tmax;
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

    Ray getRay(double u, double v, int width, int height) const {
        // Calculate camera basis vectors
        Vec3 forward = (target - position).normalize();
        Vec3 right = Vec3(forward.z, 0, -forward.x).normalize();
        Vec3 newUp = Vec3(
                forward.y * right.z - forward.z * right.y,
                forward.z * right.x - forward.x * right.z,
                forward.x * right.y - forward.y * right.x
        ).normalize();

        // Image plane dimensions
        double aspectRatio = (double)width / height;
        double viewportHeight = 2.0 * tan(fov * M_PI / 360.0);
        double viewportWidth = viewportHeight * aspectRatio;

        // Position on image plane
        double x = (u - 0.5) * viewportWidth;
        double y = (v - 0.5) * viewportHeight;

        Vec3 rayDir = forward + right * x + newUp * y;
        return Ray(position, rayDir);
    }

    Vec3 getPosition() const { return position; }
    double getRadius() const { return radius; }
    double getAngleH() const { return angleH * 180.0 / M_PI; }
    double getAngleV() const { return angleV * 180.0 / M_PI; }
};

// Function rendering the scene
void render(const Camera& camera, const Cube& cube, int width, int height) {
    vector<vector<char>> screen(height, vector<char>(width, '.'));

    // For each pixel
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            // Normalized coordinates [0, 1]
            double u = (double)x / (width - 1);
            double v = (double)(height - 1 - y) / (height - 1); // inverted Y

            Ray ray = camera.getRay(u, v, width, height);

            double t;
            if (cube.intersect(ray, t)) {
                screen[y][x] = '0';
            }
        }
    }

    // Display
    system("clear || cls"); // Clear screen (works on Linux and Windows)

    cout << "+";
    for (int i = 0; i < width; i++) cout << "-";
    cout << "+\n";

    for (int y = 0; y < height; y++) {
        cout << "|";
        for (int x = 0; x < width; x++) {
            cout << screen[y][x];
        }
        cout << "|\n";
    }

    cout << "+";
    for (int i = 0; i < width; i++) cout << "-";
    cout << "+\n";

    // Camera info
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