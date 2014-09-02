#ifndef _ARRAYIED_ALLOCATOR_H_
#define _ARRAYIED_ALLOCATOR_H_

#include <boost/static_assert.hpp>

// Visual C++ 2008以下は対応しません、GCCはC++11オプションで動かして下さい
#if (_MSC_VER <= 1500) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
	BOOST_STATIC_ASSERT(false);
#endif

#include <cstdint>

#if !defined(__INTEL_COMPILER) && !defined(__GXX_EXPERIMENTAL_CXX0X__) && (_MSC_VER < 1800)
	#include <boost/noncopyable.hpp>
#endif

namespace CheckPoint {
	template <std::size_t TTypeSize, std::size_t TNumArray>
	class ArraiedAllocator
#if !defined(__INTEL_COMPILER) && !defined(__GXX_EXPERIMENTAL_CXX0X__) && (_MSC_VER < 1800)
		: private boost::noncopyable
#endif
	{
#if defined(__INTEL_COMPILER) || defined(__GXX_EXPERIMENTAL_CXX0X__) || (_MSC_VER >= 1800)
		ArraiedAllocator(const ArraiedAllocator &) = delete;
		ArraiedAllocator & operator=(const ArraiedAllocator &) = delete;
#endif

#if defined(__INTEL_COMPILER) || defined(__GXX_EXPERIMENTAL_CXX0X__)
		static constexpr std::int32_t MAX_SIZE = TNumArray;
#else
		static const std::int32_t MAX_SIZE = TNumArray;
#endif
		// サイズは絶対０より大きくなくちゃダメ
		BOOST_STATIC_ASSERT(TNumArray > 0);

		struct Item {
			union {
				char value_[TTypeSize];
				struct Item * next_;
			};
		};

		static Item items_[MAX_SIZE];
		static Item * first_;
		static ArraiedAllocator allocator_;
		ArraiedAllocator();

	public:
		static void * Alloc() {
			Item * ret = first_;
			first_ = ret->next_;
			return reinterpret_cast<void *>(ret);
		}
		static void Free(void * item) {
			Item * rev = reinterpret_cast<Item *>(item);
			rev->next_ = first_;
			first_ = rev;
		}
		static ArraiedAllocator& GetAllocator() { return allocator_; }
		static std::int32_t Max() { return MAX_SIZE; }
	};

	template <std::size_t TTypeSize, std::size_t TNumArray>
	typename ArraiedAllocator<TTypeSize,TNumArray>::Item
		ArraiedAllocator<TTypeSize,TNumArray>::items_[ArraiedAllocator<TTypeSize,TNumArray>::MAX_SIZE];

	template <std::size_t TTypeSize, std::size_t TNumArray>
	typename ArraiedAllocator<TTypeSize, TNumArray>::Item* ArraiedAllocator<TTypeSize,TNumArray>::first_;

	template <std::size_t TTypeSize, std::size_t TNumArray>
	ArraiedAllocator<TTypeSize,TNumArray> ArraiedAllocator<TTypeSize,TNumArray>::allocator_;

	template <std::size_t TTypeSize, std::size_t TNumArray>
	inline ArraiedAllocator<TTypeSize,TNumArray>::ArraiedAllocator() {
		first_ = &items_[0];
		for (std::int32_t i = 0; i < TNumArray; i++) {
			items_[i].next_ = &items_[i + 1];
		}
		items_[TNumArray - 1].next_ = NULL;
	}
}

#endif // _ARRAYIED_ALLOCATOR_H_
