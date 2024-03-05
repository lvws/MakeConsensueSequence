#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include <algorithm>
#include <sstream>
#include <iostream>
#include "consencese.h"
#include "BamReader.h"

using std::string;

string ReverseComplement(const string &str)
{
    string result;
    result.resize(str.size());
    size_t n = str.size() - 1;
    for (const auto &base : str)
    {
        switch (base)
        {
        case 'A':
            result[n] = 'T';
            break;
        case 'T':
            result[n] = 'A';
            break;
        case 'C':
            result[n] = 'G';
            break;
        case 'G':
            result[n] = 'C';
            break;
        default:
            result[n] = 'N';
            break;
        }
        n--;
    }
    return result;
}

bool UmiMismatch1Check(const string &umi1, const string &umi2)
{
    if (umi1.size() != umi2.size())
        return false;

    // misMatch 没有影响umi 顺序
    int mis = 0;
    for (auto i = 0; i < umi1.size(); i++)
    {
        if (umi1[i] != umi2[i])
            mis++;
        if (mis > 1)
            break;
    }
    if (mis == 1)
        return true;
    else
        mis = 0;

    // misMatch 可能影响了umi 顺序
    int n = umi1.find('-');
    string umi1_a = umi1.substr(0, n);
    string umi1_b = umi1.substr(n + 1, umi1.size());
    string umi1_new = umi1_b + "-" + umi1_a;

    for (auto i = 0; i < umi1.size(); i++)
    {
        if (umi1_new[i] != umi2[i])
            mis++;
        if (mis > 1)
            return false;
    }

    return true;
}

std::tuple<std::vector<singleM>, int> UmiCluster(const std::vector<reader *> &vreads)
{
    float Threshold = 4; // 合并的比例,超过倍数的合并
    int n_umiF = 0;      // 记录错误umi 整合到簇的数量
    std::vector<singleM> v_singms;
    // 同UMI reads 成簇
    std::unordered_map<string, std::pair<std::vector<reader *>, std::vector<reader *>>> u_n2vreads_map;
    // UMI 成簇 reads 数量
    std::unordered_map<string, int> u_n2counts_map;
    int n_mapQ = 0;
    // 读取readers 记录umi 相同的，成簇
    for (auto &rd : vreads)
    {
        // std::cout << "singleM: " << rd->name << "\tvreads size: " << vreads.size() << std::endl;
        if (rd->r1Q < 15 && rd->r2Q < 15)
        {
            n_mapQ++;
            continue;
        }
        rd->umi = rd->name.substr(0, rd->name.find('|'));
        if (rd->flag1 != 0)
        {
            rd->read1 = getSeqFu8(rd->r1array, rd->r1len, rd->name);
            rd->r1q = getQualFu8(rd->q1array, rd->r1len);
        }
        else
        {
            n_mapQ++;
            continue;
        }
        if (rd->flag2 != 0)
        {
            rd->read2 = getSeqFu8(rd->r2array, rd->r2len, rd->name);
            rd->r2q = getQualFu8(rd->q2array, rd->r2len);
        }
        else
        {
            n_mapQ++;
            continue;
        }

        // std::cout << rd->read1 << std::endl;
        if (rd->flag1 & 16) // r1 reverse
        {
            rd->read1 = ReverseComplement(rd->read1);
            std::reverse(rd->r1q.begin(), rd->r1q.end());
            rd->r1_reverse = true;
        }
        if (rd->flag2 & 16) // r2 reverse
        {
            rd->read2 = ReverseComplement(rd->read2);
            std::reverse(rd->r2q.begin(), rd->r2q.end());
        }

        if (rd->r1_reverse)
            u_n2vreads_map[rd->umi].first.emplace_back(rd);
        else
            u_n2vreads_map[rd->umi].second.emplace_back(rd);
        u_n2counts_map[rd->umi] += 1;
    }

    if (n_mapQ == vreads.size())
    {
        // std::cout << "n_mapQ Failed: " << vreads[0]->name << "\t" << n_mapQ << std::endl;
        return {v_singms, n_umiF};
    }

    // 各类umi 按照 reads 数进行排序， 从多到少， 好在合并umi 时优先使用多的umi;
    std::vector<std::pair<string, int>> v_names(u_n2counts_map.begin(), u_n2counts_map.end());
    std::sort(v_names.begin(), v_names.end(), [](const std::pair<string, int> &p1, const std::pair<string, int> &p2)
              { return p1.second > p2.second; });
    for (size_t i = 0; i < v_names.size() - 1; i++)
    {
        if (v_names[i].second == 0)
            continue;
        for (size_t j = i + 1; j < v_names.size(); j++)
        {
            if (v_names[j].second == 0)
                continue;
            if (UmiMismatch1Check(v_names[i].first, v_names[j].first) && v_names[i].second > v_names[j].second * Threshold)
            {
                u_n2vreads_map[v_names[i].first].first.insert(
                    u_n2vreads_map[v_names[i].first].first.end(),
                    u_n2vreads_map[v_names[j].first].first.begin(),
                    u_n2vreads_map[v_names[j].first].first.end());
                u_n2vreads_map[v_names[i].first].second.insert(
                    u_n2vreads_map[v_names[i].first].second.end(),
                    u_n2vreads_map[v_names[j].first].second.begin(),
                    u_n2vreads_map[v_names[j].first].second.end());
                n_umiF += v_names[j].second;
                std::cout << "FIX: " << v_names[i].first << "\t" << v_names[j].first << "\n";
                v_names[j].second = 0;
            }
        }
    }

    // 构建单分子结构体
    for (auto umi : v_names)
    {
        if (umi.second != 0)
        {
            // std::cout << "Main UMI: " << umi.first << "\t";
            // std::cout << "First\t";
            // for (auto n : u_n2vreads_map[umi.first].first)
            // {
            //     std::cout << n->name << "\t";
            // }
            // std::cout << "Second\t";
            // for (auto n : u_n2vreads_map[umi.first].second)
            // {
            //     std::cout << n->name << "\t";
            // }
            // std::cout << "\n";
            v_singms.emplace_back(umi.first,
                                  u_n2vreads_map[umi.first].first,
                                  u_n2vreads_map[umi.first].second);
        }
        // else
        // {
        //     std::cout << "Fix UMI:" << umi.first << "\t";
        //     std::cout << "First\t";
        //     for (auto n : u_n2vreads_map[umi.first].first)
        //     {
        //         std::cout << n->name << "\t";
        //     }
        //     std::cout << "Second\t";
        //     for (auto n : u_n2vreads_map[umi.first].second)
        //     {
        //         std::cout << n->name << "\t";
        //     }
        //     std::cout << "\n";
        // }
    }
    return {v_singms, n_umiF};
}

/*
计算 String中 I 、D 字符出现的次数, 用于后续scs 构建筛选
*/
std::string CountID(const std::string &cigar)
{
    int i = 0, d = 0;
    for (const char &c : cigar)
    {
        switch (c)
        {
        case 'I':
            i++;
            break;
        case 'D':
            d++;
            break;
        default:
            break;
        }
    }
    return std::to_string(i) + "I" + std::to_string(d) + "D";
}

consencese MakeScs(const std::vector<reader *> &vreads, int QualThrethold, int BaseNumThrethold, float BaseThrethold)
{
    consencese n_con;
    if (vreads.size() == 0)
    {
        n_con.name = "";
        return n_con;
    }
    std::unordered_map<std::string, int> u_idcount_r1;
    std::unordered_map<std::string, int> u_idcount_r2;
    std::unordered_map<std::string, int> u_idcout2num;
    std::unordered_map<std::string, std::vector<size_t>> u_cigar2num_r1; // key: NIND(Cigar) Value: index in vector
    std::unordered_map<std::string, std::vector<size_t>> u_cigar2num_r2;
    std::cout << "Make css: ";
    for (size_t i = 0; i < vreads.size(); i++)
    {
        std::cout << vreads[i]->name << "\t";
        auto r1_ct = CountID(vreads[i]->r1cigar);
        auto r2_ct = CountID(vreads[i]->r2cigar);
        u_idcount_r1[r1_ct]++;
        u_idcount_r2[r2_ct]++;
        u_cigar2num_r1[r1_ct].push_back(i);
        u_cigar2num_r2[r2_ct].push_back(i);
    }
    std::cout << "\n";
    // 获取最主要类型的cigar。
    std::vector<std::pair<std::string, int>> r1_vpair(u_idcount_r1.begin(), u_idcount_r1.end());
    std::vector<std::pair<std::string, int>> r2_vpair(u_idcount_r2.begin(), u_idcount_r2.end());
    std::sort(r1_vpair.begin(), r1_vpair.end(), [](std::pair<std::string, int> &a, std::pair<std::string, int> &b)
              { return a.second > b.second; });
    std::sort(r2_vpair.begin(), r2_vpair.end(), [](std::pair<std::string, int> &a, std::pair<std::string, int> &b)
              { return a.second > b.second; });


    int n; // 记录 reads 数量
    char t_read[200];
    char t_qual[200];
    auto &indexes = u_cigar2num_r1[r1_vpair[0].first];
    for (size_t i = 0; i < 200; i++)
    {
        // std::cout << "N: " << n <<"\tI: " << i << std::endl;
        // std::cout << vreads[0]->read1 << std::endl;
        // std::cout << vreads[0]->r1q << std::endl;
        n = indexes.size();
        int use_base = 0;
        std::unordered_map<char, int> u_atcg_map;
        int qual = 0;
        for (auto j : indexes)
        {
            if (i < (vreads[j]->read1).size())
            {
                if ((vreads[j]->r1q)[i] >= QualThrethold)
                {
                    u_atcg_map[(vreads[j]->read1)[i]]++;
                    qual += (vreads[j]->r1q)[i];
                    use_base += 1;
                }
            }
            else
            {
                n--;
            }
        }
        // std::cout << "2:N: " << n <<"\tI: " << i << std::endl;
        if (n == 0) // 没有多余的碱基
        {
            t_read[i] = 0;
            t_qual[i] = 0;
            break;
        }
        // ============ 碱基比较 ============
        t_read[i] = 'N';
        t_qual[i] = '!';
        if (use_base < BaseNumThrethold)
        {
            for (const auto &c : u_atcg_map)
            {
                // std::cout << "u_atcg_map: " << "\t" << c.second << std::endl;
                // std::cout << "base is : " << c.first << std::endl;
                if (c.second == use_base)
                {
                    t_read[i] = c.first;
                    t_qual[i] = qual / use_base;
                    break;
                }
            }
        }
        else
        {
            for (const auto &c : u_atcg_map)
            {
                if (c.second >= use_base * BaseThrethold)
                {
                    t_read[i] = c.first;
                    t_qual[i] = qual / use_base;
                    break;
                }
            }
        }
        // for (const auto &c : u_atcg_map)
        // {
        //     std::cout << c.first << ": " << c.second << " Threthold: " << use_base * BaseThrethold
        //               << "\tuse_base: " << use_base << "\tBaseThrethold: " << BaseThrethold << std::endl;
        // }
        // std::cout << "use: " << t_read[i] << std::endl;
        // std::cout << "3:N: " << n <<"\tI: " << i << std::endl;
        //
    }
    // std::cout << "SCS-Con: "  << std::endl;
    n_con.r1 = t_read;
    n_con.q1 = t_qual;
    // std::cout << "SCS-Con-A: "  << std::endl;
    /*
    ============= 处理R2 ==============
    */
    char t_read2[200];
    char t_qual2[200];
    indexes = u_cigar2num_r2[r2_vpair[0].first];
    for (size_t i = 0; i < 200; i++)
    {
        n = indexes.size();
        int use_base = 0;
        std::unordered_map<char, int> u_atcg_map;
        int qual = 0;
        for (auto j : indexes)
        {
            if (i < (vreads[j]->read2).size())
            {
                if ((vreads[j]->r2q)[i] >= QualThrethold)
                {
                    u_atcg_map[(vreads[j]->read2)[i]]++;
                    qual += (vreads[j]->r2q)[i];
                    use_base++;
                }
            }
            else
            {
                n--;
            }
        }

        if (n == 0) // 没有多余的碱基
        {
            t_read2[i] = 0;
            t_qual2[i] = 0;
            break;
        }
        t_read2[i] = 'N';
        t_qual2[i] = '!';
        if (use_base < BaseNumThrethold)
        {
            for (const auto &c : u_atcg_map)
            {
                if (c.second == use_base)
                {
                    t_read2[i] = c.first;
                    t_qual2[i] = qual / use_base;
                    break;
                }
            }
        }
        else
        {
            for (const auto &c : u_atcg_map)
            {
                if (c.second >= use_base * BaseThrethold)
                {
                    t_read2[i] = c.first;
                    t_qual2[i] = qual / use_base;
                    break;
                }
            }
        }
    }
    // std::cout << "SCS-Con: "  << std::endl;
    n_con.r2 = t_read2;
    n_con.q2 = t_qual2;
    n_con.name = vreads[0]->name;
    n_con.id1 = r1_vpair[0].first;
    n_con.id2 = r2_vpair[0].first;
    // std::cout << "SCS-Con-B: "  << std::endl;

    return n_con;
}

int IdNumCount(const string &s)
{
    std::istringstream ss(s);
    int i, d;
    char _i, _d;
    ss >> i >> _i >> d >> _d;
    return i + d;
}

consencese MakeDcs(const consencese &con1, const consencese &con2)
{
    consencese n_con;
    if (con1.name == "")
        return con2;
    if (con2.name == "")
        return con1;
    n_con.name = con1.name;
    n_con.s1 = con1.s1;
    n_con.s2 = con2.s2;
    int con1_id, con2_id;

    // 处理R1
    if (con1.id1 != con2.id2) // 某一条序列可能有错误插入缺失
    {
        con1_id = IdNumCount(con1.id1);
        con2_id = IdNumCount(con2.id2);
        if (con1_id < con2_id)
        {
            n_con.r1 = con1.r1;
            n_con.q1 = con1.q1;
        }
        else
        {
            n_con.r1 = con2.r2;
            n_con.q1 = con2.q2;
        }
    }
    else
    {
        char t_read[200];
        char t_qual[200];

        for (size_t i = 0; i < 200; i++)
        {
            int n_use = 0;
            if (i < con1.r1.size())
            {
                t_read[i] = con1.r1[i];
                t_qual[i] = con1.q1[i];
                n_use++;
            }
            if (i < con2.r2.size())
            {
                if (n_use == 0)
                {
                    t_read[i] = con2.r2[i];
                    t_qual[i] = con2.q2[i];
                }
                else
                {
                    if (t_read[i] != con2.r2[i])
                    {
                        t_read[i] = 'N';
                        t_qual[i] = '!';
                    }
                    else
                    {
                        t_qual[i] = (t_qual[i] + con2.q2[i]) / 2;
                    }
                }
                n_use++;
            }

            if (n_use == 0)
            {
                t_read[i] = 0;
                t_qual[i] = 0;
                break;
            }
        }
        n_con.r1 = t_read;
        n_con.q1 = t_qual;
    }

    /*
    ================= R2 ===============
    */
    if (con1.id2 != con2.id1) // 某一条序列可能有错误插入缺失
    {
        con1_id = IdNumCount(con1.id2);
        con2_id = IdNumCount(con2.id1);
        if (con1_id < con2_id)
        {
            n_con.r2 = con1.r2;
            n_con.q2 = con1.q2;
        }
        else
        {
            n_con.r2 = con2.r1;
            n_con.q2 = con2.q1;
        }
    }
    else
    {
        char t_read2[200];
        char t_qual2[200];
        for (size_t i = 0; i < 200; i++)
        {
            int n_use = 0;
            if (i < con1.r2.size())
            {
                t_read2[i] = con1.r2[i];
                t_qual2[i] = con1.q2[i];
                n_use++;
            }
            if (i < con2.r1.size())
            {
                if (n_use == 0)
                {
                    t_read2[i] = con2.r1[i];
                    t_qual2[i] = con2.q1[i];
                }
                else
                {
                    if (t_read2[i] != con2.r1[i])
                    {
                        t_read2[i] = 'N';
                        t_qual2[i] = '!';
                    }
                    else
                    {
                        t_qual2[i] = (t_qual2[i] + con2.q1[i]) / 2;
                    }
                }
                n_use++;
            }

            if (n_use == 0)
            {
                t_read2[i] = 0;
                t_qual2[i] = 0;
                break;
            }
        }
        n_con.r2 = t_read2;
        n_con.q2 = t_qual2;
    }

    return n_con;
}
