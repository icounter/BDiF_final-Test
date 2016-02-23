///create by Chao Zhu
#include <iostream>
#include "mpi.h"
#include <sstream>
#include <vector>
#include "boost/date_time/posix_time/posix_time.hpp"
#include <numeric>
#include <math.h>
#include <string.h>
#include <unordered_map>
#include <functional>
#include <queue>
#include "main.h"
using namespace std;
using namespace boost::posix_time;
using namespace boost::gregorian;

ptime time_t_epoch(date(1970, 1, 1));
struct record{
    date  date1;
    ptime time;
    float price;
    int volume;
};
vector<string> splitline(char* buff){
    string a=string(buff);
    vector<string> result;
    int j=0;
    for(int i=0;i<a.size();++i){
        if(a[i]==',') {
            result.push_back(a.substr(j,i-j));
            j=i+1;
        }
    }
    if(a[a.size()-1]!=',') result.push_back(a.substr(j,a.size()-j));
//    for(int i=0;i<result.size();i++){
//        cout<<result[i]<<endl;
//    }
    return result;
}
pair<float,float> calculate_mean_sd(vector<float> resultSet){
    double sum = std::accumulate(std::begin(resultSet), std::end(resultSet), 0.0);
    float mean =  sum / resultSet.size(); //均值

    float accum  = 0.0;
    std::for_each (std::begin(resultSet), std::end(resultSet), [&](const float d) {
        accum  += (d-mean)*(d-mean);
    });

    float stdev = sqrt(accum/(resultSet.size()-1)); //方差
    return make_pair(mean,stdev);
}

vector<float> stat(vector<record> buff){
    vector<float> price;
    vector<float> volume;
    for(int i=0;i<buff.size();++i){
        price.push_back(log(buff[i].price));
        volume.push_back(log(float(buff[i].volume)));
    }
    for(int i=0;i<volume.size();++i){
        if(volume[i]<0) cout<<buff[i].volume<<endl;
    }
    vector<float> result;
    pair<float,float> pri,volu;
    pri=calculate_mean_sd(price);
    volu=calculate_mean_sd(volume);
    float price_mean=pri.first;
    float price_sd=pri.second;
    float vol_mean=volu.first;
    float vol_sd=volu.second;
    result.push_back(exp(price_mean-4*price_sd));
    result.push_back(exp(price_mean+4*price_sd));
    result.push_back(vol_mean-3*vol_sd);
    result.push_back(vol_mean+3*vol_sd);
    return result;
}
struct record convert_string_to_record(vector<string> line){
    struct record record1;
    date d(from_undelimited_string(line[0].substr(0,8)));
    string a=" ";
    string ts(to_simple_string(d)+" "+line[0].substr(9,15));
    ptime t(time_from_string(ts));
    float price(stof(line[1]));
    int volume(stoi(line[2]));
    record1.date1=d;
    record1.time=t;
    record1.price=price;
    record1.volume=volume;
    return record1;
}
bool Compare(pair<string,int> a1,pair<string,int> a2){
    return a1.second<a2.second;
}


void process(char *buff,vector<record> &valid,vector<record> &noise){
    char *pch;
    pch=strtok(buff,"\n");
    int i=0;
    vector<record> temp_vector;
    while(pch!=NULL){
        vector<string> line=splitline(pch);
        if(line.size()==3&&line[0].size()==24){//this is potential right answer
            struct record temp=convert_string_to_record(line);
            if(temp.price<=0||temp.volume<=0||temp.price>100000||temp.volume>10000000){
                noise.push_back(temp);
            }
            else {
                temp_vector.push_back(temp);
                i++;
                if (i > 10000) {
                    i = i - 10000;
                    vector<float> min_max_number = stat(temp_vector);
                   // cout<<min_max_number[0]<<endl;
                    priority_queue<pair<string, int>, vector<pair<string, int>>, function<bool(pair<string, int>,
                                                                                           pair<string, int>)>> count_date(
                            Compare);
                    unordered_map<string,int> map;
                    for(int i=0;i<temp_vector.size();i++){
                        unordered_map<string,int>::iterator got = map.find (to_simple_string(temp_vector[i].date1));
                        if(got==map.end()){
                            map.insert(make_pair(to_simple_string(temp_vector[i].date1),1));
                        }
                        else{
                            got->second++;
                        }
                    }
                    for (auto it = map.begin(); it != map.end(); ++it){
                        count_date.push(make_pair(it->first,it->second));
                    }
                    unordered_map<string, int> datemap;
                    datemap.insert({count_date.top().first,count_date.top().second});
                    count_date.pop();
                    if (count_date.size() >= 1 && count_date.top().second > 100) {
                        datemap.insert({count_date.top().first,count_date.top().second});
                    }//if the second largest number is >100 is reasonable
                    for (int i = 0; i < temp_vector.size(); i++) {//&&temp_vector[i].volume > exp(min_max_number[2]) &&temp_vector[i].volume < exp(min_max_number[3])
                        unordered_map<string,int>::iterator got = datemap.find(to_simple_string(temp_vector[i].date1));
                        if (temp_vector[i].price > min_max_number[0] && temp_vector[i].price < min_max_number[1] &&
                            got!=datemap.end()) {
                            valid.push_back(temp_vector[i]);
                        }
                        else {
                            noise.push_back(temp_vector[i]);
                        }

                    }
                    temp_vector.clear();
                }
            }
        }
        pch = strtok (NULL, "\n");
    }
    if(i!=0&&temp_vector.size()!=0){
        vector<float> min_max_number=stat(temp_vector);
        cout<<min_max_number[0]<<"  "<<min_max_number[1]<<min_max_number[2]<<min_max_number[3]<<endl;
        priority_queue<pair<string,int>,vector<pair<string,int>>,function<bool(pair<string,int>,pair<string,int>)>> count_date(Compare);
        unordered_map<string,int> map;
        for(int i=0;i<temp_vector.size();i++){
            unordered_map<string,int>::iterator got = map.find (to_simple_string(temp_vector[i].date1));
            if(got==map.end()){
                map.insert(make_pair(to_simple_string(temp_vector[i].date1),1));
            }
            else{
                got->second++;
            }
        }
        for (auto it = map.begin(); it != map.end(); ++it){
            count_date.push(make_pair(it->first,it->second));
        }
        unordered_map<string, int> datemap;
        datemap.insert({count_date.top().first,count_date.top().second});
        count_date.pop();
        if (count_date.size() >= 1 && count_date.top().second > 100) {
            datemap.insert({count_date.top().first,count_date.top().second});
        }//if the second largest number is >100 is reasonable
        for(int i=0;i<temp_vector.size();i++){//&&temp_vector[i].volume>exp(min_max_number[2])&&temp_vector[i].volume<exp(min_max_number[3])
            unordered_map<string,int>::iterator got = datemap.find(to_simple_string(temp_vector[i].date1));
            if(temp_vector[i].price>min_max_number[0]&&temp_vector[i].price<min_max_number[1]&&
              got!=datemap.end()){
                valid.push_back(temp_vector[i]);
            }
            else{
                noise.push_back(temp_vector[i]);
            }
        }
        temp_vector.clear();
    }
}
string convert_ptime_to_string(ptime time){
    string a=to_iso_extended_string(time);
    string b;
    b=a.substr(0,4)+a.substr(5,2)+a.substr(8,2)+':'+a.substr(11,8)+'.'+a.substr(20,6);
    return b;

}
int64_t GetCurrentStamp64()
{
    boost::posix_time::ptime epoch(boost::gregorian::date(1970, boost::gregorian::Jan, 1));
    boost::posix_time::time_duration time_from_epoch =
              boost::posix_time::microsec_clock::universal_time() - epoch;
    return time_from_epoch.total_microseconds();

}
/////////////normality

void normal_test(vector<record> valid,float &param,long long ndim){
vector<float> vec;
    for(int i=0;i<valid.size();i++){
        vec.push_back(log(valid[i].price));
    }
    Normality norm(vec);
    param=norm.get_JB_factor();
    ndim=vec.size();

}









int main(int argc, char** argv) {
    int rank, size;
    long long read_cnt;
    long long filesize;
    MPI_File fh,fhv,fhn,fhl;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    char* FILENAME=argv[1];
    int64_t t1=GetCurrentStamp64();
    int ret= MPI_File_open(MPI_COMM_WORLD,FILENAME,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);
    MPI_File_get_size(fh,&filesize);
    long long local_size=filesize/size;
    read_cnt=(rank==size-1)?((filesize-(size-1)*local_size)):local_size;
    char *buff=(char *)(malloc(sizeof(char)*(read_cnt+1)));
    if(buff==NULL) {cout<<"error!"; return 0;}
  MPI_File_read_at(fh, rank*local_size, buff, read_cnt, MPI_BYTE, &status);
    int64_t t2=GetCurrentStamp64();

    vector<record> valid;
    vector<record> noise;
   process(buff,valid,noise);
    int64_t t3=GetCurrentStamp64();
   char *valid_out=new char[2*read_cnt];
//   for(int i=0;i<valid.size();i++){
//       cout<<to_simple_string(valid[i].time)<<","<<to_string(valid[i].price)<<","<<to_string(valid[i].volume)<<endl;
////    }
//    for(int i=0;i<noise.size();i++){
//        cout<<to_simple_string(noise[i].time)<<","<<to_string(noise[i].price)<<","<<to_string(noise[i].volume)<<endl;
//    }
   //cout<<"noise.size"<<noise.size()<<endl;
char* p;
    p=valid_out;
    for(int i=0;i<valid.size();i++){
        string a=convert_ptime_to_string(valid[i].time)+','+to_string(valid[i].price)+','+to_string(valid[i].volume);
        strcpy(p,a.c_str());
        p+=a.length();
        *p='\n';
        p++;
    }
    char *noise_out=new char[read_cnt];
    char *q;
    q=noise_out;
    for(int i=0;i<noise.size();i++){
        string b;
        b=convert_ptime_to_string(noise[i].time)+','+to_string(noise[i].price)+','+to_string(noise[i].volume);
        strcpy(q,b.c_str());
        q+=b.length();
        *q='\n';
        q++;
   }
//    strcpy(noise_out,b.c_str());
    MPI_File_open(MPI_COMM_WORLD,"valid_data.txt",MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fhv);
    MPI_File_set_size(fhv,0);
    MPI_File_open(MPI_COMM_WORLD,"noise_data.txt",MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fhn);
    MPI_File_set_size(fhn,0);
   //cout<<"reter"<<retr<<" "<<endl;
  //  MPI_File_close(&fhv);
    long valid_allsz[size];
    long write_valid_size=strlen(valid_out);
    long noise_allsz[size];
    long write_noise_size=strlen(noise_out);
    MPI_Allgather(&write_valid_size,1,MPI_LONG,valid_allsz,1,MPI_LONG,MPI_COMM_WORLD);
    MPI_Allgather(&write_noise_size,1,MPI_LONG,noise_allsz,1,MPI_LONG,MPI_COMM_WORLD);
    for(int i=1;i<size;i++){
        valid_allsz[i]+=valid_allsz[i-1];
    }
    for(int i=1;i<size;i++){
        noise_allsz[i]+=noise_allsz[i-1];
    }
    cout<<write_noise_size<<endl;
    //cout<<valid_allsz[rank]<<endl;
    MPI_File_write_at_all(fhv,valid_allsz[rank]-write_valid_size,valid_out,write_valid_size,MPI_BYTE,&status);
    MPI_File_write_at_all(fhn,noise_allsz[rank]-write_noise_size,noise_out,write_noise_size,MPI_BYTE,&status);
    //splitfile(buff);//plitfile buff
//for(int i=0;i<read_cnt;++i){
//    cout<<buff[i];
//}

//cout<<"||||||||||"<<endl;
//    long write_size=read_cnt/2;
//    int allsz[size];

    free(noise_out);

    free(buff);
    MPI_File_close(&fhv);
    MPI_File_close(&fhn);
    MPI_File_close(&fh);
    int64_t t4=GetCurrentStamp64();
    stringstream ss;
    ss<<"rank:"<<rank<<" "<<"time for read data "<<(t2-t1)<<endl;
    ss<<"rank:"<<rank<<" "<<"time for process data "<<t3-t2<<endl;
    ss<<"rank:"<<rank<<" "<<"time for write data "<<t4-t3<<endl;
    string sss=ss.str();
    char * buffer=new char [sss.length()];
    strcpy(buffer,sss.c_str());
    MPI_File_open(MPI_COMM_WORLD,"logfile.txt",MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fhl);
    MPI_File_set_size(fhl,0);
    int write_time=sss.length();
    int time_allsz[size];
    MPI_Allgather(&write_time,1,MPI_INT,time_allsz,1,MPI_INT,MPI_COMM_WORLD);
    for(int i=1;i<size;i++){
        time_allsz[i]+=time_allsz[i-1];
    }
    MPI_File_write_at_all(fhl,time_allsz[rank]-write_time,buffer,write_time,MPI_BYTE,&status);
    free(buffer);
    MPI_File_close(&fhl);
    float JB_factor;
    long long ndim;
    normal_test(valid,JB_factor,ndim);
    stringstream norm;
    norm<<"rank: "<<rank<<"JB stats is"<<JB_factor<<"in how many data points"<<ndim<<endl;
    string norm1=norm.str();
    char* buffer1=new char [norm1.length()];
    strcpy(buffer1,norm1.c_str());
    MPI_File_open(MPI_COMM_WORLD,"normtest.txt",MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fhl);
    MPI_File_set_size(fhl,0);
    int norm_time=norm1.length();
    int norm_allsz[size];
    MPI_Allgather(&norm_time,1,MPI_INT,norm_allsz,1,MPI_INT,MPI_COMM_WORLD);
    for(int i=1;i<size;i++){
        norm_allsz[i]+=norm_allsz[i-1];
    }
    MPI_File_write_at_all(fhl,norm_allsz[rank]-norm_time,buffer1,norm_time,MPI_BYTE,&status);
    free(buffer1);
    MPI_File_close(&fhl);
    free(valid_out);
    MPI_Finalize();
    return 0;
















}

