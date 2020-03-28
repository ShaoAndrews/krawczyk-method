#pragma once

template <class T>
class Node
{
public:
    T region;
    Node<T> *child1;
    Node<T> *child2;
    Node<T> *child3;
    Node<T> *child4;

public:
    Node(T region_);
    ~Node();
};

template <class T>
Node<T>::Node(T region_) : region(region_)
{
    child1 = nullptr;
    child2 = nullptr;
    child3 = nullptr;
    child4 = nullptr;
}
template <class T>
Node<T>::~Node()
{
    delete child1;
    delete child2;
    delete child3;
    delete child4;
    child1 = nullptr;
    child2 = nullptr;
    child3 = nullptr;
    child4 = nullptr;
}
