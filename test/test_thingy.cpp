#include <iostream>
#include <memory>
#include <set>

struct Tester {
  int t = 0;
};

using shared = std::shared_ptr<Tester>;
using weak = std::weak_ptr<Tester>;

bool operator<(weak lhs, weak rhs) {
  return lhs.lock() < rhs.lock();
}

int main() {
  unsigned i = 32;
  std::cout << i << std::endl;
  std::cout << (i << 2) << std::endl;
  std::cout << (i >> 2) << std::endl;

  return 0;
}
