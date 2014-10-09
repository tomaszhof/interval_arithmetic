/*
 * THashMap.cpp
 *
 *  Created on: Dec 29, 2012
 *      Author: tomaszhof
 */

#include "THashMap.h"

using namespace interval_arithmetic;

namespace interval_arithmetic {

template<typename T>
THashMap<T>::THashMap() {
	// TODO Auto-generated constructor stub

}

template<typename T>
THashMap<T>::~THashMap() {
	// TODO Auto-generated destructor stub
}

template<typename T>
Interval<T> THashMap<T>::FromMap(int ij) {
	Interval<T> tmp = { 0, 0 };
	if (map.find(ij) != map.end()) {
		tmp = map[ij];
	}

	return tmp;
}

template<typename T>
void THashMap<T>::ToMap(int ij, Interval<T> v) {
	if (map.find(ij) != map.end()) {
		map[ij] = v;
	} else
		map.insert(std::map<int, Interval<T> >::value_type(ij, v));
}

template<typename T>
void THashMap<T>::Clear() {
	map.clear();
}

template<typename T>
void THashMap<T>::Erase(int ij) {
	map.erase(ij);
}
}
