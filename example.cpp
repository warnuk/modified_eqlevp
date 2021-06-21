#include <iostream>
using namespace std;


// Experimenting with classes and scope
class example {
    public:

    int x = 0;
    int y = 0;

    example(int z) {
        x += z;
        y += z;
    }

    void test(int a) {
        x *= a;
        y *= a;
    }
};

int main() {
    example x1(10);
    cout << x1.x << endl;

    x1.test(3);

    cout << x1.x << endl;

    return 0;
}