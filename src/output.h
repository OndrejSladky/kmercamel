#pragma once

#include <iostream>
#include <string>
#include <vector>

/// Print the fasta file header.
void WriteName(const int k) {
    std::cout << ">superstring ";
    std::cout << "k=" << k << std::endl;
}

/// Mask the superstring and print it to the standard output.
void WriteSuperstring(std::string superstring, std::vector<bool> mask) {
    std::string ret;
    ret.reserve(superstring.length());
    for (size_t i = 0; i < superstring.length(); ++i) {
        ret +=  mask[i] ? (char)superstring[i] : (char)std::tolower(superstring[i]);
    }
    std::cout << ret;
}

