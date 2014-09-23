/*
 * THashMap.h
 *
 *  Created on: Dec 29, 2012
 *      Author: tomaszhof
 */

#ifndef THASHMAP_H_
#define THASHMAP_H_

#include "Interval.h"
#include <map>

using namespace interval_arithmetic;

namespace interval_arithmetic
{

template<typename T>
class THashMap
{
private:
	std::map<int, Interval<T>> map;
public:
	THashMap();
	virtual ~THashMap();
	void ToMap(int ij, Interval<T> v);
	Interval<T> FromMap(int ij);
	void Clear();
	void Erase(int ij);
};

}

#endif /* THASHMAP_H_ */
