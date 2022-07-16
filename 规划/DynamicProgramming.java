package leetcode;

import org.apache.groovy.json.internal.Value;

import java.lang.reflect.Array;
import java.util.Arrays;

public class DynamicProgramming {

    /**
     * 未分类/分硬币/状态确定
     * <p>
     * leetcode # 322<p>
     * 给你一个整数数组 coins ，表示不同面额的硬币；以及一个整数 amount ，表示总金额。<p>
     * 计算并返回可以凑成总金额所需的 最少的硬币个数 。如果没有任何一种硬币组合能组成总金额，返回-1 。<p>
     * 你可以认为每种硬币的数量是无限的。<p>
     * 来源：力扣（LeetCode）<p>
     * 链接：https://leetcode-cn.com/problems/coin-change<p>
     *
     * @param coin:硬币包含的面值
     * @return :可以凑成总金额所需的 最少的硬币个数 。如果没有任何一种硬币组合能组成总金额，返回-1 。
     */
    public static int coin(int[] coin, int aim) {
        int[] dp = new int[aim + 1];
        /**
         * dp[i]表示凑齐 i 需要的最少硬币数
         * dp[i] = MAX_value 表示无法凑出金额 i
         * i<0 时dp[i] = MAX_value
         */
        dp[0] = 0;
        for (int i = 1; i <= aim; i++) {
            /**
             * Integer.MAX_VALUE + 1 = - Integer.MAX_VALUE
             * Integer.MAX_VALUE-1 防止超出整形范围
             */
            dp[i] = Integer.MAX_VALUE - 1;
            int min = Integer.MAX_VALUE - 1; //防止超出整形范围
            for (int j = 0; j < coin.length; j++) {
                if (i - coin[j] < 0) {
                    min = Math.min(Integer.MAX_VALUE, min);
                } else {
                    min = Math.min(min, dp[i - coin[j]]);
                }
            }
            dp[i] = min + 1;
        }
        return dp[aim] < Integer.MAX_VALUE ? dp[aim] : -1;
    }


    /**
     * 坐标型/跳跃/状态<p>
     * leetcode #55<p>
     * 给定一个非负整数数组 nums ，你最初位于数组的 第一个下标 。<p>
     * 数组中的每个元素代表你在该位置可以跳跃的最大长度。<p>
     * 判断你是否能够到达最后一个下标。<p>
     *
     * @param nums:非负整数数组
     * @return : true or false
     */
    public boolean canJump(int[] nums) {
        boolean[] dp = new boolean[nums.length + 1];
        /**
         * dp[i] 表示能否可以跳到 位置i
         */
        dp[0] = true;
        dp[1] = true;
        for (int i = 2; i <= nums.length; i++) {
            dp[i] = false;
            for (int j = 1; j < i; j++) {
                if (dp[j] && i - j <= nums[j - 1]) {
                    dp[i] = true;
                    break;
                }
            }
        }
        return dp[nums.length];
    }

    /**
     * 坐标型/跳跃/最值 TODT<p>
     * leetcode #45<p>
     * 给你一个非负整数数组nums ，你最初位于数组的第一个位置。<p>
     * 数组中的每个元素代表你在该位置可以跳跃的最大长度。<p>
     * 你的目标是使用最少的跳跃次数到达数组的最后一个位置。<p>
     * 假设你总是可以到达数组的最后一个位置。<p>
     *
     * @param nums：非负整数数组
     * @return ：最少的跳跃次数
     */
    public int jump(int[] nums) {
        /**
         * dp[i] 表示跳到 i的最小跳跃次数
         */
        int[] dp = new int[nums.length];

        dp[0] = 0;
        for (int i = 1; i < dp.length; i++) {
            dp[i] = nums.length + 1;
        }

        for (int i = 0; i < nums.length; i++) {
            for (int j = 1; j <= nums[i]; j++) {
                if (i + j >= nums.length) {
                    return dp[dp.length - 1];
                }
                dp[i + j] = Math.min(dp[i + j], dp[i] + 1);// 从前到后
            }
        }

        return dp[dp.length - 1];
    }

    /**
     * 坐标型/网格移动2/最值<p>
     * leetcode #64<p>
     * 给定一个包含非负整数的 m x n 网格 grid ，请找出一条从左上角到右下角的路径，使得路径上的数字总和为最小。<p>
     * 说明：每次只能向下或者向右移动一步。<p>
     *
     * @param grid:
     * @return :
     */
    public int minPathSum(int[][] grid) {
        int[][] dp = new int[2][grid[0].length];
        /**
         * dp[i][k] 表示走到 grid[i][k] 的路径上的数字最小总和
         * dp[i][k] = min( dp[i-1][k] + grid[i][k], dp[i][k-1] + grid[i][k] )
         */
        dp[0][0] = grid[0][0];

        /**
         * 数组轮换
         */
        int now = 1, old = 0;
        int up, left;
        for (int i = 0; i < grid.length; i++) {
            old = now;
            now = (now + 1) % 2;// old now 在0，1间轮换

            for (int j = 0; j < grid[0].length; j++) {
                if (i == 0 && j == 0)
                    continue;

                if (i < 1) {
                    up = Integer.MAX_VALUE;
                } else {
                    up = dp[old][j] + grid[i][j];
                }

                if (j < 1) {
                    left = Integer.MAX_VALUE;
                } else {
                    left = dp[now][j - 1] + grid[i][j];
                }
                dp[now][j] = Math.min(up, left);
            }
        }
        return dp[now][grid[0].length - 1];

    }


    /**
     * 序列型/盗窃房屋1/最值<p>
     * leetcode #198<p>
     * 相邻的房屋装有相互连通的防盗系统，<p>
     * 如果两间相邻的房屋在同一晚上被小偷闯入，系统会自动报警。<p>
     * 给定一个代表每个房屋存放金额的非负整数数组，<p>
     * 计算不触动警报装置的情况下 ，一夜之内能够偷窃到的最高金额。<p>
     *
     * @param nums：每个房屋存放金额的非负整数数组
     * @return ：偷窃到的最高金额
     */
    public int rob(int[] nums) {
        int[][] dp = new int[nums.length + 1][2];
        /**
         * dp[i][k]表示前 i栋获得的最大值，下标为 i-1 房屋的状态为 k
         */
        dp[0][0] = 0;
        dp[0][1] = 0;
        for (int i = 1; i <= nums.length; i++) {
            for (int k = 0; k <= 1; k++) {
                if (k == 0) {
                    dp[i][k] = Math.max(dp[i - 1][0], dp[i - 1][1]);
                } else {
                    //k == 1
                    dp[i][k] = dp[i - 1][0] + nums[i - 1];
                }
            }
        }
        return Math.max(dp[nums.length][0], dp[nums.length][1]);
    }

    /**
     * 序列型/paint house/最值<p>
     * lintcode # 516<p>
     * 这里有n个房子在一列直线上，现在我们需要给房屋染色，<p>
     * 共有k种颜色。每个房屋染不同的颜色费用也不同，你希望每两个相邻的房屋颜色不同<p>
     * 费用通过一个nxk 的矩阵给出，比如cost[0][0]表示房屋0染颜色0的费用，<p>
     * cost[1][2]表示房屋1染颜色2的费用。找到油漆所有房子的最低成本。<p>
     *
     * @param costs:
     * @return :
     */
    public static int paintHouseMinCost(int[][] costs) {
        if (costs.length == 0) {
            return 0;
        }

        int[][] dp = new int[costs.length + 1][costs[0].length];
        /**
         * dp[i][k] 表示前 i栋的最小花费，第 i栋的颜色为 k
         *
         */
        System.arraycopy(costs[0], 0, dp[1], 0, costs[0].length);

        for (int i = 2; i <= costs.length; i++) {

            int min = Integer.MAX_VALUE, min_2 = Integer.MAX_VALUE;
            int min_index = 0, min_index2 = 0;
            for (int t = 0; t < dp[i - 1].length; t++) {
                if (dp[i - 1][t] < min_2) {
                    min_2 = dp[i - 1][t];
                    min_index2 = t;
                    if (min_2 < min) {
                        int temp = min;
                        min = min_2;
                        min_2 = temp;
                        min_index2 = min_index;
                        min_index = t;
                    }
                }
            }//找到最小值和次小值的坐标

            for (int k = 0; k < costs[0].length; k++) {
                if (k == min_index) {
                    dp[i][k] = dp[i - 1][min_index2] + costs[i - 1][k];
                } else {
                    dp[i][k] = dp[i - 1][min_index] + costs[i - 1][k];
                }
            }
        }

        Arrays.sort(dp[costs.length]);
        return dp[costs.length][0];
    }


    /**
     * 序列型/盗窃房屋2/最值<p>
     * leetcode #213<p>
     * 相邻的房屋装有相互连通的防盗系统，<p>
     * 如果两间相邻的房屋在同一晚上被小偷闯入，系统会自动报警。<p>
     * 这个地方所有的房屋都 围成一圈 ，这意味着第一个房屋和最后一个房屋是紧挨着的。<p>
     * 给定一个代表每个房屋存放金额的非负整数数组，<p>
     * 计算不触动警报装置的情况下 ，一夜之内能够偷窃到的最高金额。<p>
     *
     * @param nums：每个房屋存放金额的非负整数数组
     * @return ：偷窃到的最高金额
     */
    public int rob2(int[] nums) {
        if (nums.length <= 1) {
            return nums[0];
        }

        int[][][] dp = new int[nums.length + 1][2][2];
        /**
         * dp[i][j][k]表示前 i栋获得的最大值，下标 为 i-1 房屋的状态为 j,第 一 栋的状态为 k
         */
        dp[1][1][1] = nums[0];
        for (int i = 2; i <= nums.length; i++) {
            for (int j = 0; j <= 1; j++) {
                for (int k = 0; k <= 1; k++) {
                    if (j == 0) {
                        dp[i][j][k] = Math.max(dp[i - 1][0][k], dp[i - 1][1][k]);
                    } else {
                        //j == 1
                        dp[i][j][k] = dp[i - 1][0][k] + nums[i - 1];
                    }
                }

            }
        }
        return Math.max(dp[nums.length][0][0], Math.max(dp[nums.length][1][0], dp[nums.length][0][1]));
    }

    /**
     * 序列型+位运算<p>
     * 给你一个整数 n ，对于 0 <= i <= n 中的每个 i ，<p>
     * 计算其二进制表示中 1 的个数 ，<p>
     * 返回一个长度为 n + 1 的数组 ans 作为答案。<p>
     *
     * @param n：
     * @return ：
     */
    public int[] countBits(int n) {
//        int[] dp = new int[n + 1];
//        dp[0] = 0;
//        for (int i = 1; i <= n; i++) {
//            dp[i] = Integer.bitCount(i);
//        }
//        return dp;
        int[] dp = new int[n + 1];
        dp[0] = 0;
        if (n == 0) {
            return dp;
        }
        for (int i = 1; i <= n; i++) {
            dp[i] = dp[i >> 1] + i % 2;
        }
        return dp;
    }

    /**
     * 序列型/子序列问题/最长严格递增子序列的长度<p>
     * leetcode # 300<p>
     *
     * 一个整数数组 nums ，找到其中最长严格递增子序列的长度。<p>
     * 子序列是由数组派生而来的序列，删除（或不删除）数组中的元素而不改变其余元素的顺序。<p>
     * 例如，[3,6,2,7] 是数组 [0,3,1,6,2,2,7] 的子序列。<p>
     *
     * @param nums：
     * @return ：<p>
     *
     *
     * 考虑数组上常用的两种思路：<p>
     *
     * 每次减少一半：如果每次将问题规模减少一半，原问题有[10,9,2,5]，和[3,7,101,18]，<p>
     * 两个子问题的最优解分别为 [2,5] 和 [3,7,101]，但是找不到好的组合方式将两个子问题最优解组合为原问题最优解 [2,5,7,101]。<p>
     *
     * 每次减少一个：记 f(n) 为以第 n 个数结尾的最长子序列，每次减少一个，<p>
     * 将原问题分为 f(n-1), f(n-2), ..., f(1)，共 n - 1 个子问题。n−1=7 个子问题以及答案如下：<p>
     * [10, 9, 2, 5, 3, 7, 101] -> [2, 5, 7, 101]<p>
     * [10, 9, 2, 5, 3, 7] -> [2, 5, 7]<p>
     * [10, 9, 2, 5, 3] -> [2, 3]<p>
     * [10, 9, 2, 5] -> [2, 5]<p>
     * [10, 9, 2] -> [2]<p>
     * [10, 9] -> [9]<p>
     * [10] -> [10]<p>
     *
     * 已经有 7 个子问题的最优解之后，可以发现一种组合方式得到原问题的最优解：f(6)f(6) 的结果 [2,5,7],<p>
     * 7 < 18，同时长度也是 f(1)~f(7) 中，结尾小于 18 的结果中最长的。<p>
     * f(7) 虽然长度为 4 比 f(6) 长，但结尾是不小于 18 的，无法组合成原问题的解。<p>
     *
     * 以上组合方式可以写成一个式子，即状态转移方程
     */
    public static int lengthOfLIS(int[] nums) {
        int[] dp = new int[nums.length];
        /**
         * dp[n] 表示以第 n+1 个数结尾的最长子序列
         * dp[n] = max(dp[i]) + 1  , i < n , nums[i] < nums[n]
         */
        dp[0] = 1;
        for (int i = 1; i < nums.length; i++) {
            int max = 0;
            for (int k = 0; k < i; k++) {
                if (nums[i] > nums[k]) {
                    max = Math.max(max, dp[k]);
                }
            }
            dp[i] = Math.max(max, dp[i]) + 1;
        }
        Arrays.sort(dp);
        return dp[nums.length - 1];
    }

    /**
     * 序列型/子序列/最大和的连续子数组<p>
     * leetcode #<p>
     * 给你一个整数数组 nums ，请你找出一个具有最大和的连续子数组（子数组最少包含一个元素），返回其最大和。<p>
     * 子数组 是数组中的一个连续部分。<p>
     *
     * @param nums:
     * @return :
     */
    public int maxSubArray(int[] nums) {
        int[] dp = new int[nums.length];
        /**
         * dp[n] 表示以 nums[n] 为最后一个元素的连续子序列的最大和值
         * dp[n] = max(dp[n-1], 0) + nums[n]
         */
        dp[0] = nums[0];
        for (int i = 1; i < nums.length; i++) {
            dp[i] = Math.max(dp[i - 1], 0) + nums[i];
        }
        Arrays.sort(dp);
        return dp[nums.length - 1];
    }


    public static void main(String[] args) {
        int[] nums = {10, 9, 2, 1};
        System.out.println(lengthOfLIS(nums));
    }


}
