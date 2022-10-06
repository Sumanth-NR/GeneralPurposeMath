#include <bits/stdc++.h>
using namespace std;

class Solution {
public:
	vector<int> searchRange(vector<int> &nums, int target) {
		int l = lower_bound(nums.begin(), nums.end(), target) - nums.begin();
		if (l == nums.size()) return { -1, -1 };
		int r = upper_bound(nums.begin(), nums.end(), target + 1) - nums.begin() - 1;
		return { l, r };
	}
};

int main() {
	int t; cin >> t;
	Solution s;
	vector<int> ans = s.searchRange({5, 7, 7, 8, 8, 10}, 8);
	cout << ans[0] << ' ' << ans[1] << '\n';
	return 0;
}
