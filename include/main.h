#include <string>
#include <chrono>
#include "BamReader.h"
#include "consencese.h"

using std::string;

/*
计时器
*/
class Timer
{
private:
    std::chrono::_V2::system_clock::time_point now;
public:
    Timer(): now(std::chrono::high_resolution_clock::now()) {}

    ~Timer()
    {
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<float> duration = end - now;
        float ms = duration.count() * 1000.0f;
        std::cout << "Timer: " << ms << " ms" << std::endl;
    }

};

/*
容器，用于快速存储、删除数据，需要预先设定大小。
*/
template <typename T>
struct MyVector
{
private:
    T* data;
    size_t start,end,len,num;
public:
    MyVector(unsigned n)
        : start(0),end(0),len(n),num(0)
    {
        data = new T[n];
    }

    ~MyVector()
    {
        delete[] data;
    }

    bool push_back(T elem)
    {
        if (num == len)
            return false;
        data[end] = elem;
        end ++ ;
        num ++ ;
        if (end == len) 
            end = 0;
        return true;
    }

    T pop()
    {
        T tmp;
        if (num == 0)
            return tmp;
        tmp = data[start];
        start ++ ;
        num -- ;
        if (start == len) 
            start = 0;
        return tmp;
    }

    // 检查头部元素
    T& head()
    {
        return data[start];
    }

    // 检查尾部元素
    T& tail()
    {
        return data[end];
    }

    // 返回容器的容量
    size_t size()
    {
        return len;
    }

    // 返回容器中目前元素数量
    size_t elems()
    {
        return num;
    }

    void status()
    {
        std::cout << "MyVector 容量：" << len  << "; ";
        std::cout << "存储元素：" << num << "; " ;
        std::cout << "起始指针位置：" << start << ";" ;
        std::cout << "终止指针位置：" << end << "\n";
    }
};

/*
从输入的BAM， aln 信息中，构建reader, 会修改 reader 
*/
void ParseSamRead(const bam1_t *aln, reader* rd);

/*
读取BAM文件，构建reader 、consences
*/
void ReadBAM(const string& bam_name, std::unordered_map<string,reader*>& u_name2reader_map,
std::unordered_map<string,std::vector<reader*>> u_pos2readers_map,
std::vector<std::string> v_splitreads);


