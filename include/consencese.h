#include <string>
#include <vector>
#include <iostream>
#include <memory>

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
    int flag1; // r1 flag 
    int flag2; // r2 flag
    bool r1_reverse = false; // 记录r1 是否为反向map 到基因组，好对 cc 进行DCS构建

    uint8_t *r1array,*q1array,*r2array,*q2array; // 序列、质量的 数组
    int r1len,r2len; // 序列长度

    // ~reader(){
    //     delete[] r1array;
    //     delete[] r2array;
    //     delete[] q1array;
    //     delete[] q2array;
    // }

};

/*
consencese 数据结构，方便后续写出fastq
*/
struct consencese
{
    string r1,r2;
    string q1,q2;
    string id1,id2; // del and insert ; 记录del 和 ins 的状态.
    string name;
    int s1 = 0;
    int s2 = 0;
    /* data */
};

/*
UMI 单分子聚类数据结构
*/
struct singleM
{
    string umi;
    std::vector<reader*> vr1s ; // r1 反向比对的单分子链
    std::vector<reader*> vr2s ; // r1 正向比对的单分子链

    singleM(const string& umi, const std::vector<reader*>& vr1s, const std::vector<reader*>& vr2s )
        :umi(umi),vr1s(vr1s),vr2s(vr2s){}
};

// 序列的反向互补
string ReverseComplement(const string& str);

/*
检查UMI之间的MisMatch，
由于开头碱基的变化会导致，UMI排序变化，所以正反比较两次
*/
bool UmiMismatch1Check(const string& umi1, const string& umi2);

/*
对同一个起始终止区域的所有reads 按照umi 进行聚类；
misMatch 1 ,多UMI 以最主要的umi 为准（ > 90% ?）
*/
std::vector<singleM> UmiCluster(const std::vector<reader*>& vreads);

/*
构建单链一致性序列，
*/
consencese MakeScs(const std::vector<reader*>& vreads,int QualThrethold,int BaseNumThrethold,float BaseThrethold);

/*
构建双链一致性序列，
*/
consencese MakeDcs(const consencese& con1, const consencese& con2);
