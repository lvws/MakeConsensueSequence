#include <string>
#include <unordered_map>
#include <vector>

using std::string;

/*
read pair 的存储格式，用于记录BAM中的read pair ，以进行后续consenses构建 
*/
struct reader
{
    string name;
    string read1; // 所有的read序列还原到正常顺序
    string read2;
    string r1q;
    string r2q;
    string r1cigar; // cigar值用来后续构建consences 插入缺失的参考，防止错位导致大量的不匹配
    string r2cigar;
    string umi; // ACATGACT-CTGCT
    string umi_tag; // ab
    bool r1_reverse; // 记录r1 是否为反向map 到基因组，好对 cc 进行DCS构建
};
// 序列的反向互补
string ReverseComplement(const string& str);



