#include "greedy.cpp"

std::vector<KMer> ReadData() {
    int count;
    std::cin >> count;
    std::vector<KMer> out(count);
    for(int i = 0; i < count; ++i) {
        std::cin >> out[i].value;
    }
    return out;
}

void WriteResult(KMerSet result) {
    std::cout << result.superstring << std::endl;
    for (auto x : result.mask) {
        std::cout << x;
    }
    std::cout << std::endl;
}

int main() {
    auto data = ReadData();
    auto result = Greedy(data);
    WriteResult(result);
    return 0;
}
