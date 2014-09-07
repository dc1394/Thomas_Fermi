#ifndef _MKL_ALLOCATOR_H_
#define _MKL_ALLOCATOR_H_

#ifdef _MSC_VER
	#pragma once
#endif

#include <new>
#include <limits>
#include <cstddef>
#include <cstdint>
#include <mkl.h>

namespace thomasfermi {
	template <typename T>
	class mkl_allocator
	{
	public:
		// typedefs
		typedef std::size_t size_type;
		typedef std::ptrdiff_t difference_type;
		typedef T * pointer;
		typedef const T * const_pointer;
		typedef T & reference;
		typedef const T & const_reference;
		typedef T value_type;
			
		// convert an allocator<T> to allocator<U>
		template <class U> 
		struct rebind { typedef mkl_allocator<U> other; };

		// constructors
		mkl_allocator() throw() {}
		mkl_allocator(mkl_allocator const &) throw() {}
		template <typename U>
		mkl_allocator(mkl_allocator<U> const &) throw() {}

		// destructor
		~mkl_allocator() throw() {}

		// メモリを割り当てる
		pointer allocate(size_type size, const_pointer hint = 0) {
			void * const p = mkl_malloc(size * sizeof(T), 64);
			if (!p)
				throw std::bad_alloc();

			return reinterpret_cast<pointer>(p);
		}

		// 割当て済みの領域を初期化する
		void construct(pointer p, const T & val)
		{ new (reinterpret_cast<void *>(p)) T(val);	}

		// メモリを解放する
		void deallocate(pointer p, size_type n)
		{ mkl_free(reinterpret_cast<void *>(p)); }

		// 初期化済みの領域を削除する
		void destroy(pointer p)
		{ p->~T(); }

		// アドレスを返す
		pointer address(reference value) const
		{ return &value; }
		const_pointer address(const_reference value) const
		{ return &value; }
		
		// 割当てることができる最大の要素数を返す
		size_type max_size() const throw()
		{ return (std::numeric_limits<std::size_t>::max()) / sizeof(T); }
	};

	template <typename T, typename U>
	inline bool operator ==(const mkl_allocator<T>&, const mkl_allocator<U>)
	{ return true; }

	template <typename T, typename U>
	inline bool operator !=(const mkl_allocator<T>&, const mkl_allocator<U>)
	{ return false; }
}

#endif	// _MKL_ALLOCATOR_H_
