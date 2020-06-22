/*
 * Matrix2D.h
 *
 *  Created on: 22 cze 2020
 *      Author: tomhof
 */

#ifndef SRC_MATRIX2D_H_
#define SRC_MATRIX2D_H_

// A dynamic size matrix using std::vector for storage.

//--------------------------------------------- Machinery:
#include <algorithm>        // std::copy
#include <assert.h>         // assert
#include <initializer_list> // std::initializer_list
#include <vector>           // std::vector
#include <stddef.h>         // ptrdiff_t

using Size = ptrdiff_t;
using std::initializer_list;
using std::vector;

template< class Item >
class Matrix2D
{
private:
	vector<Item>    items_;
	Size            n_cols_;

	auto index_for( Size const x, Size const y ) const
		-> Size
	{ return y*n_cols_ + x; }

public:
	auto n_rows() const -> Size { return items_.size()/n_cols_; }
	auto n_cols() const -> Size { return n_cols_; }

	auto item( Size const x, Size const y )
		-> Item&
	{ return items_[index_for(x, y)]; }

	auto item( Size const x, Size const y ) const
		-> Item const&
	{ return items_[index_for(x, y)]; }

	Matrix2D(): n_cols_( 0 ) {}

	Matrix2D( Size const n_cols, Size const n_rows )
		: items_( n_cols*n_rows )
		, n_cols_( n_cols )
	{}

	Matrix2D( initializer_list< initializer_list<Item> > const& values )
		: items_()
		, n_cols_( values.size() == 0? 0 : values.begin()->size() )
	{
		for( auto const& row : values )
		{
			assert( Size( row.size() ) == n_cols_ );
			items_.insert( items_.end(), row.begin(), row.end() );
		}
	}
};




#endif /* SRC_MATRIX2D_H_ */
