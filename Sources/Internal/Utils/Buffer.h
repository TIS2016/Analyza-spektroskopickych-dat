//
//  Buffer.h
//  Utils
//
//  Created by Filip Kovac on 19/10/2016.
//  Copyright © 2016 Filip Kovac. All rights reserved.
//

#include "MemoryAllocator.h"

template < class T, class Allocator = MemoryAllocator<T> > class Buffer
{
private:
	size_t	m_size;
	T*			m_pItems;
	
public:
	
	inline T* Ptr() { return m_pItems; }
	
	inline size_t Length() { return m_size; }
	
	HRESULT Allocate ( size_t length )
	{
		T *temp = m_pItems;
		if ( temp != nullptr )
		{
			Allocator::Release( temp );
		}
		
		m_pItems = Allocator::Allocate( length );
		if ( m_pItems == nullptr )
		{
			return E_OUTOFMEMORY;
		}
		m_size = length;
		
		return S_OK;
	}
	
};
