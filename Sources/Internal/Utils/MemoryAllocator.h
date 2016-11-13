//
//  MemoryAllocator.h
//  Utils
//
//  Created by Filip Kovac on 07/10/2016.
//  Copyright © 2016 Filip Kovac. All rights reserved.
//

#include <stdlib.h>

namespace Internal {
	
#pragma push_macro("new")
#pragma push_macro("delete")

#undef new
#undef delete

template < class T > class MemoryAllocatorImplementation {
private:
	struct AllocationHelper {
		T o;
		
		void* operator new ( size_t size )
		{
			return malloc( size );
		}
		
		void* operator new[] ( size_t size )
		{
			return malloc( size );
		}
		
		void operator delete ( void *ptr )
		{
			free( ptr );
		}
		
		void operator delete[] ( void *ptr )
		{
			free( ptr );
		}
	};
	
public:
	
	T* Allocate ( size_t size )
	{
		return (T*)(new AllocationHelper[size]);
	};
	
	T* AllocateItem( size_t size )
	{
		return (T*)(new AllocationHelper(size));
	}
	
	void Release( T* buffer )
	{
		delete[] (AllocationHelper *)(buffer);
	}
	
	void ReleaseItem( T* pItem )
	{
		delete (AllocationHelper *)(pItem);
	}
	
	
};

#pragma pop_macro("new")
#pragma pop_macro("delete")
	
}

template < class T > class MemoryAllocator: public Internal::MemoryAllocatorImplementation< T > {};
