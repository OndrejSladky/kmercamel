#include "greedy.cpp"
#include "kmers.cpp"

#include <iostream>
#include <string>

constexpr int k = 20;

std::string ReadData() {
    std::string result;
    std::cin >> result;
    return result;
}

void WriteResult(KMerSet result, std::vector<KMer> kMers, std::string data) {
    std::cout << result.superstring << std::endl;
    for (auto x : result.mask) {
        std::cout << x;
    }
    std::cout << std::endl;
    std::cout << "superstring length:         " << result.superstring.length() << std::endl;
    std::cout << "k-mers count:               " << kMers.size() << std::endl;
    std::cout << "length of scanned sequence: " << data.length() << std::endl;
    std::cout << "coefficient:                " << result.superstring.length() / (double)kMers.size() << std::endl;

}

int main() {
    auto data = ReadData();
    auto kMers = ConstructKMers(data, k);
    auto result = Greedy(kMers);
    WriteResult(result, kMers, data);
    return 0;
}
