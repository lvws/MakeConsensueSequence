#include <iostream>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <tuple>
#include <fstream>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <vector>
#include <queue>
#include <memory>

#include "main.h"

using std::cout, std::cerr, std::string;

// 线程安全的队列
template<typename T>
class ThreadSafeQueue {
private:
    std::queue<T> queue_;
    mutable std::mutex mutex_;
    std::condition_variable cond_;
    int num = 0;

public:
    ThreadSafeQueue() {}

    void push(T value) {
        std::lock_guard<std::mutex> lock(mutex_);
        queue_.push(value);
        cond_.notify_one();
        num++;
    }

    bool try_pop(T& value) {
        std::lock_guard<std::mutex> lock(mutex_);
        if (queue_.empty()) {
            return false;
        }
        value = queue_.front();
        queue_.pop();
        num--;
        return true;
    }

    void wait_and_pop(T& value) {
        std::unique_lock<std::mutex> lock(mutex_);
        cond_.wait(lock, [this]() { return !queue_.empty(); });
        value = queue_.front();
        queue_.pop();
        num--;
    }

    bool empty() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return queue_.empty();
    }

    int elems() {
        return num;
    }
};

// 非处理测试
void FreeReads(std::vector<reader*> data)
{
    for (auto rd : data){
        delete rd;
    }
}

// 单线程测试
void MakeCs1(std::vector<reader*> data, const std::string& file_name) {
    int QualThrethold = 10 + 33;
    int BaseNumThrethold = 3;
    float BaseThrethold = 0.75;
    std::ofstream output1(file_name+"_R1.fastq", std::ios::out); // 以追加模式打开文件
    std::ofstream output2(file_name+"_R2.fastq", std::ios::out); 

    // std::cout << "开始取数据：" << std::endl;
    // 从队列中取出数据
    
    //std::cout << "数据量：" << data.size() << std::endl;
    // 检查是否为特殊的终止值
    if (data.size() == 0) {
        return;
    }
    // std::cout << "Start: MakeCs" << std::endl;
    auto v_umi_clusters = UmiCluster(data);
    // std::cout << "umi_clusters is OK!"<< std::endl;
    std::vector<consencese> v_out_cons;
    for (auto& sm : v_umi_clusters)
    {
        //std::cout << "SCS1: " << sm.vr1s.size() << "\tSCS2: " << sm.vr2s.size()  << std::endl;      
        auto con1 = MakeScs(sm.vr1s,QualThrethold,BaseNumThrethold,BaseThrethold);
        con1.s1 = sm.vr1s.size();
        std::cout << "make scs1 ok!"<< std::endl;
        auto con2 = MakeScs(sm.vr2s,QualThrethold,BaseNumThrethold,BaseThrethold);
        con2.s2 = sm.vr2s.size();
        std::cout << "make scs2 ok!"<< std::endl;
        v_out_cons.emplace_back(MakeDcs(con1,con2));
        std::cout << "make Dcs ok!"<< std::endl;

        // 清空已经用完的reads 数据
        for (auto i: sm.vr1s)
        {
            std::cout << "预备删除vr1s：" << i->name << std::endl;
        }
        for (auto i: sm.vr2s)
        {
            std::cout << "预备删除vr2s：" << i->name << std::endl;
        }
        for (auto i: sm.vr1s)
        {
            delete i;
        }
        for (auto i: sm.vr2s)
        {
            delete i;
        }
    }

    // 用互斥锁保护文件写入
    {
        for (const auto& con : v_out_cons)
        {
            output1 << con.name << "|" << con.s1 << "|" << con.s2 << "\n";
            output1 << con.r1 << "\n";
            output1 << "+\n";
            output1 << con.q1 << "\n";

            output2 << con.name << "|" << con.s1 << "|" << con.s2 << "\n";
            output2 << con.r2 << "\n";
            output2 << "+\n";
            output2 << con.q2 << "\n";
        }
        
    }
    std::cout << "Finish MakeCs1!" << std::endl;

}


// 工作线程的函数
void MakeCs(ThreadSafeQueue<std::vector<reader*>>& data_queue, std::mutex& file_mutex, 
const std::string& file_name, int threadNum,
int QualThrethold,int BaseNumThrethold,float BaseThrethold) 
{
    std::ofstream output1(file_name+"_R1.fastq"+std::to_string(threadNum), std::ios::out); 
    std::ofstream output2(file_name+"_R2.fastq"+std::to_string(threadNum), std::ios::out); 
    std::vector<reader*> data;
    while (true) {
        //std::cout << "开始取数据：" << std::endl;
        // 从队列中取出数据
        data_queue.wait_and_pop(data);
        //std::cout << "数据量：" << data.size() << std::endl;
        // 检查是否为特殊的终止值
        if (data.size() == 0) {
            break;
        }
        //std::cout << "Start: MakeCs" << std::endl;
        auto v_umi_clusters = UmiCluster(data);
        //std::cout << "umi_clusters is OK!"<< std::endl;
        std::vector<consencese> v_out_cons;
        for (auto& sm : v_umi_clusters)
        {
            //std::cout << "SCS1: " << sm.vr1s.size() << "\tSCS2: " << sm.vr2s.size()  << std::endl;
            auto con1 = MakeScs(sm.vr1s,QualThrethold,BaseNumThrethold,BaseThrethold);
            con1.s1 = sm.vr1s.size();
            //std::cout << "make scs1!"<< std::endl;
            auto con2 = MakeScs(sm.vr2s,QualThrethold,BaseNumThrethold,BaseThrethold);
            con2.s2 = sm.vr2s.size();
            //std::cout << "make scs2!"<< std::endl;
            v_out_cons.emplace_back(MakeDcs(con1,con2));
            //std::cout << "make Dcs!"<< std::endl;

            // 清空已经用完的reads 数据
            for (auto i: sm.vr1s)
            {
                delete i;
            }
            for (auto i: sm.vr2s)
            {
                delete i;
            }
        }

        // 用互斥锁保护文件写入
        {
            // std::lock_guard<std::mutex> lock(file_mutex);
            for (const auto& con : v_out_cons)
            {
                //std::cout << con.name << std::endl;
                output1 << "@" << con.name << "|" << con.s1 << "|" << con.s2 << "\n";
                output1 << con.r1 << "\n";
                output1 << "+\n";
                output1 << con.q1 << "\n";

                output2 << "@" << con.name << "|" << con.s1 << "|" << con.s2 << "\n";
                output2 << con.r2 << "\n";
                output2 << "+\n";
                output2 << con.q2 << "\n";
            }
        }
    }
}

void ParseSamRead(const bam1_t *aln, reader* rd)
{
    int flag = getFlag(aln);
    rd->name = getName(aln);

    if (flag & 64) // fist in pair
    {
        rd->r1cigar = getCigar(aln);
        rd->read1 = getSeq(aln);
        rd->r1q = getQual(aln);
        rd->flag1 = flag;
    } else {
        rd->r2cigar = getCigar(aln);
        rd->read2 = getSeq(aln);
        rd->r2q = getQual(aln);
        rd->flag2 = flag;
    }
}

void ReadBAM(const string& bam_name,ThreadSafeQueue<std::vector<reader*>>& data_queue)
{
    int MaxInsertSize = 10000;
    /*
    方便查询同一个起始位点reads集合
    */
    MyVector<std::tuple<string,string,int>> v_readspoition(MaxInsertSize); 
    /*
    根据read name 找到对应的read
    */
    std::unordered_map<string,reader*> u_name2reader_map;

    /*
    read pair pos-insertsize map， 将同起点、终点的read pairs 记录在一个字典里
    key = chr:pos-insert
    value: point to reader
    */
    std::unordered_map<string,std::vector<reader*>> u_pos2readers_map;
    /*
    map 在不同区域的read pairs 单独处理
    */
    std::vector<std::string> v_splitreads;
    BamReader myBam(bam_name);
    std::unordered_map<std::string,bool> u_already_recode2;

    while (myBam.next())
    {
        bool exist_tag = false; // 记录read 是否存在
        bool exist_pos_tag = false; // 记录pos-key 是否在 u_pos2readers_map 中
        bool recode_tag = false; // 记录read 是否已经记录到 u_pos2readers_map 中
        string pos2readers_key; // 记录比对位置；
        string name = getName(myBam.aln);
        int flag = getFlag(myBam.aln);
        reader* newRd;
    
        if (flag & 2304) // not primary alignment or supplementary alignment
            continue;
        // 该read 曾经记录过
        exist_tag = (u_name2reader_map.find(name) != u_name2reader_map.end());
        // 该read  已经记录到 u_pos2readers_map 中
        recode_tag = (u_already_recode2.find(name) != u_already_recode2.end());
        if (recode_tag)
            u_already_recode2.erase(name); // 这个对内存影响很大。
        // std::cout << "Start Treat: " << name << "\tExists: " << exist_tag << std::endl;
        if (exist_tag)
        {
            newRd = u_name2reader_map[name];
            ParseSamRead(myBam.aln, newRd);
            u_name2reader_map.erase(name);
        } else {
            newRd = new reader; // 在堆分配数据，好跨作用域使用，用完要记得delete!!!
            //reader* newRd = std::make_shared<reader>();
            u_name2reader_map[name] = newRd;
            ParseSamRead(myBam.aln, newRd);
        }

        // std::cout << "ParseSamRead OK! " << name << std::endl;
        string chrom = getChrom(myBam.aln,myBam.bam_header);
        int pos = getPos(myBam.aln);
        int pos2 = getNpos(myBam.aln);
        int insize = getIsize(myBam.aln);

        //  不记录insert size为负的reads , 避免重复。
        if (insize < 0)
            continue;

        if (insize > 0) // 正常比对的read pairs , 取出现在靠前的reads 
        {
            pos2readers_key = chrom + ":" + std::to_string(pos) + ":" + std::to_string(insize);
            exist_pos_tag = (u_pos2readers_map.find(pos2readers_key) != u_pos2readers_map.end());
            // std::cout << "Key: " << pos2readers_key << "\tExists: " << exist_pos_tag << std::endl; 
            if (!recode_tag) // 该read 第一出现，则加入u_pos2readers_map 中
            {
                u_pos2readers_map[pos2readers_key].emplace_back(newRd);
                u_already_recode2[name] = true;
            }
            // 准备就绪的reads 进行cs 构建
            while (true && v_readspoition.elems() > 0 )
            {
                auto& [h_pos_key,h_chrom,h_pos] = v_readspoition.head();
                if ((chrom == h_chrom && pos > h_pos) || chrom != h_chrom)
                {
                    // std::cout << "Deal with: " << h_pos_key << std::endl;  
     
                    data_queue.push(u_pos2readers_map[h_pos_key]);
                    // FreeReads(u_pos2readers_map[h_pos_key]);
                    u_pos2readers_map.erase(h_pos_key);
                    // if (u_already_recode.find(h_pos_key) != u_already_recode.end())
                    //     std::cout << h_pos_key << "\t dup" << std::endl;
                    // MakeCs1(u_pos2readers_map[h_pos_key],"output");
                    v_readspoition.pop();
                } else {
                    break;
                }
            }
            
        }
        else if(insize == 0) // 存在split reads，记录在一起，后面一并处理。
        {
            string chrom2 = getNchrom(myBam.aln,myBam.bam_header);
            pos2readers_key = chrom + ":" + std::to_string(pos) + "-" + chrom2 + ":" + std::to_string(pos2);
            exist_pos_tag = (u_pos2readers_map.find(pos2readers_key) != u_pos2readers_map.end());
            if (!recode_tag) // 该read 第一出现，则加入u_pos2readers_map 中
            {
                u_pos2readers_map[pos2readers_key].emplace_back(u_name2reader_map[name]) ;
                u_already_recode2[name] = true;
            } 
            
        }

        if (!exist_pos_tag) 
            {
                // 有些reads 的insert size 特别大
                if (insize < MaxInsertSize*0.8 && insize > 0)
                {
                    // std::cout << "Push: " << pos2readers_key << std::endl;
                    v_readspoition.push_back({pos2readers_key,chrom,std::max(pos+insize,pos2)});
        
                } else  { // 由于 insert size < 0 已经被过滤这里肯定是大于等于0；
                    v_splitreads.emplace_back(pos2readers_key);
                }
            }
    }
    std::cout << "Finish good reads!" << std::endl;
    while (v_readspoition.elems() > 0)
    {
        auto[h_pos_key,h_chrom,h_pos] = v_readspoition.pop();
        data_queue.push(u_pos2readers_map[h_pos_key]);
        // MakeCs1(u_pos2readers_map[h_pos_key],"output");
        // FreeReads(u_pos2readers_map[h_pos_key]);
        u_pos2readers_map.erase(h_pos_key);
    }

    std::cout << "Finish v_readspoition!" << std::endl;
    std::cout << "Deal v_splitreads size: " << v_splitreads.size() << std::endl;
    for (auto& key : v_splitreads)
    {
        data_queue.push(u_pos2readers_map[key]);
        // MakeCs1(u_pos2readers_map[key],"output");
        // for (auto& i :u_pos2readers_map[key] )
        //     std::cout << i->name << "\n";
        // FreeReads(u_pos2readers_map[key]);
        u_pos2readers_map.erase(key);
    }
        
}

int main() {
    const int num_threads = 10;
    ThreadSafeQueue<std::vector<reader*>> data_queue;
    std::mutex file_mutex;
    std::vector<std::thread> pool;
    std::string bam_name = "test.bam";
    const std::string file_name = "output";
    int QualThrethold = 10 + 33;
    int BaseNumThrethold = 3;
    float BaseThrethold = 0.75;
    // std::ofstream output1(file_name+"_R1.fastq", std::ios::out); 
    // std::ofstream output2(file_name+"_R2.fastq", std::ios::out); 
    
    // 启动线程池
    for (int i = 0; i < num_threads; ++i) {
        pool.emplace_back(MakeCs, std::ref(data_queue), std::ref(file_mutex),file_name,i,
        QualThrethold,BaseNumThrethold,BaseThrethold);
    }
    
    // 主线程生成数据
    // for (int i = 0; i < 100; ++i) {
    //     data_queue.push(i);
    // }
    ReadBAM(bam_name,data_queue);

    // 发送终止信号给工作线程
    for (int i = 0; i < num_threads; ++i) {
        data_queue.push({}); // 使用-1表示没有更多的数据了
    }

    // 等待所有线程完成
    for (auto& thread : pool) {
        thread.join();
    }

    // 程序结束
    //std::cout << "All data processed and written to file." << std::endl;
    return 0;
}
