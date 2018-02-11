//
//  set.h
//  example
//
//  Created by 张鑫禹 on 2017/10/29.
//
//

#ifndef set_h
#define set_h
#include <iostream>
//#define MAX 2147483647
using namespace std;

template<typename T>
class set{
private:
    T *item;
    long number;
    long ntotal;
public:
    set(long total){
        this->number = 0;
        ntotal = total;
        item = new T[total];
        memset(this->item, -1, sizeof(T)*ntotal);
    }
    /*~set(){
        delete[] item;
    }*/
    int isExist(const T item);
    bool addItem(const T item);
    bool removeItem(const T item);
    set operator+ (set set2);
    set operator* (set set2);
    void display();
    int scale();
    T visitItem(int index);
    bool isEmpty();
    T end();
};

template<typename T>
int set<T>::isExist(const T item){
    for(int i = 0;i < this->number;i++){
        if(this->item[i] == item)
            return i;
    }
    return -1;
}

template<typename T>
bool set<T>::addItem(const T item){
    if(isExist(item)>=0 || this->number>=ntotal-1){
        cerr<<"the adding itme is out of boundary";
        return false;
    }
    this->item[this->number++] = item;
    return true;
}

template<typename T>
bool set<T>::removeItem(const T item){
    int pos = isExist(item);
    if(pos<0){
        cout<<"the itme is not found";
        return false;
    }
    this->number--;
    for(int i = pos;i < this->number;i++){
        this->item[i] = this->item[i+1];
    }
    return true;
}

template<typename T>
set<T> set<T>::operator+ (set<T> set2){  //union
    set<T> result(2*this->ntotal);
    for(int i = 0;i < this->number;i++){
        result.addItem(this->item[i]);
    }
    for(int i = 0;i < set2.number;i++){
        if(result.isExist(set2.item[i]) == -1)
            result.addItem(set2.item[i]);
    }
    return result;
}

template<typename T>
set<T> set<T>::operator* (set<T> set2){  //intersection
    set<T> result(this->ntotal);
    for(int i = 0;i < this->number;i++){
        if(set2.isExist(this->item[i]) >= 0){
            result.addItem(this->item[i]);
        }
    }
    return result;
}

template<typename T>
void set<T>::display(){
    for(int i = 0;i < this->number;i++){
        cout<<this->item[i]<<" ";
    }
    cout<<"\n"<<endl;
}

template<typename T>
int set<T>::scale(){
    return number;
}

template<typename T>
T set<T>::visitItem(int index){
    if(index<0 || index >= this->number){
        cerr<<"the visiting index is out of boundary";
        exit(-1);
    }
    else return this->item[index];
}

template<typename T>
bool set<T>::isEmpty(){
    return (this->number <= 0) ? true : false;
}

template<typename T>
T set<T>::end(){
    return this->item[this->number-1];
}
#endif /* set_h */
