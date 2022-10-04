#include <string>
#include <vector>
#include <iostream>
#include <queue>
#include <unordered_set>

struct KMer {
	std::string value;

	size_t length() { return value.size(); }
};

/// Find the length of longest overlap in O(length^2).
int OverlapLength(KMer first, KMer second) {
	for (int start = 0; start < first.length(); ++start) {
		bool overlaps = true;
		for (int i = 0; i < std::min(first.length() - start, second.length()); ++i) {
			if (first.value[start + i] != second.value[i]) {
				overlaps = false;
				break;
			}
		}
		if (overlaps) return first.length() - start;
	}
	return 0;
}

/// Construct the union of two k-mers in O(length).
KMer Union(KMer first, KMer second, int overlapLength) {
	return KMer{
		first.value + second.value.substr(overlapLength, second.length() - overlapLength)
	};
}

struct KMerSet {
	std::string superstring;
	std::vector<bool> mask;
	int k;
};

struct UnionEvent {
	int firstIndex;
	int secondIndex;
	int overlapLength;
};

/// Get the approximated shortest superstring of the given k-mers using the GREEDY algorithm.
/// This runs in O(n^2 o^2) where o is the size of the output (n <= o <= nk).
KMerSet Greedy(std::vector<KMer> input) {
	if (input.size() == 0) {
		throw new std::invalid_argument("input cannot be empty");
	}
	int k = input[0].length();

	auto kMers = std::vector<KMer>(input);
	// There will be in total n - 1 new k-mers.
	std::vector<bool> used (2 * kMers.size() - 1);
	auto compare = [](UnionEvent &a, UnionEvent &b) {return a.overlapLength < b.overlapLength;};
	std::priority_queue<UnionEvent, std::vector<UnionEvent>, decltype(compare)> q(compare);
	// Find the overlaps of the initial k-mers.
    for (int i = 0; i < kMers.size(); ++i) {
		for (int j = 0; j < kMers.size(); ++j) if (j != i) {
			q.push(UnionEvent{i, j, OverlapLength(kMers[i], kMers[j])});
		}
	} 

    // Unite the k-mers hungrily.
	// After n - 1 unions, all k-mers are united.
	while (kMers.size() != 2 * input.size() - 1) {
		auto next = q.top();
		q.pop();
		// If none of the individual strings have yet been merged, merge them.
		if (!used[next.firstIndex] && !used[next.secondIndex]) {
			used[next.firstIndex] = used[next.secondIndex] = true;
			KMer united = Union(kMers[next.firstIndex], kMers[next.secondIndex], next.overlapLength);
            // Get the overlap of the newly created string with the existing ones.
            for (int i = 0; i < kMers.size(); ++i) if (!used[i]) {
				q.push(UnionEvent{i, (int)kMers.size(), OverlapLength(kMers[i], united)});
				q.push(UnionEvent{(int)kMers.size(), i, OverlapLength(united, kMers[i])});
			}
			kMers.push_back(united);
		}
	}

	auto resultValue = kMers[kMers.size() - 1].value;
	KMerSet result = KMerSet {
		resultValue,
		std::vector<bool> (resultValue.size()),
		k
	};

    // Fill in the mask.
	// TODO: make this non-specific to the representation.
	std::unordered_set<std::string> inputValues;
	for (auto x : input) inputValues.insert(x.value);
	for (int i = 0; i <= resultValue.size() - k; ++i) {
		if (inputValues.count(resultValue.substr(i, k))) {
			result.mask[i] = true;
		}
	}

	return result;
}


