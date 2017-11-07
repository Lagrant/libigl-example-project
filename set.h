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

using namespace std;

class set{
private:
    int *item;
    long number;
    long ntotal;
public:
    set(long total){
        this->number = 0;
        ntotal = total;
        item = new int[total];
        memset(this->item,0,sizeof(item));
    }
    int isExist(const int item);
    bool addItem(const int item);
    bool removeItem(const int item);
    set operator+ (set set2);
    set operator* (set set2);
    void display();
};

int set::isExist(const int item){
    for(int i = 0;i < number;i++){
        if(this->item[i] == item)
            return i;
    }
    return -1;
}

bool set::addItem(const int item){
    if(isExist(item)>=0 || this->number>=ntotal-1)
        return false;
    this->item[this->number++] = item;
    return true;
}

bool set::removeItem(const int item){
    int pos = isExist(item);
    if(pos<0)
        return false;
    this->number--;
    for(int i = pos;i < this->number;i++){
        this->item[i] = this->item[i+1];
    }
    return true;
}

set set::operator+ (set set2){  //union
    set result(2*this->ntotal);
    for(int i = 0;i < this->number;i++){
        result.addItem(this->item[i]);
    }
    for(int i = 0;i < set2.number;i++){
        if(result.isExist(set2.item[i]) == -1)
            result.addItem(set2.item[i]);
    }
    return result;
}

set set::operator* (set set2){  //intersection
    set result(this->ntotal);
    for(int i = 0;i < this->number;i++){
        if(set2.isExist(this->item[i]) >= 0){
            result.addItem(this->item[i]);
        }
    }
    return result;
}

void set::display(){
    for(int i = 0;i < this->number;i++){
        cout<<this->item[i]<<" ";
    }
    cout<<"\n"<<endl;
}

#endif /* set_h */
