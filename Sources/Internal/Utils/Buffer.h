//
//  Buffer.h
//  Utils
//
//  Created by Filip Kovac on 19/10/2016.
//  Copyright © 2016 Filip Kovac. All rights reserved.
//

#include "stdafx.h"
#include "MemoryAllocator.h"

namespace DataAnalysis { namespace Utils {

	template < class T, class Allocator = MemoryAllocator<T> > class Buffer
	{
	private:
		size_t	m_size;
		T*		m_pItems;
	
	public:
	
		Buffer () {
			m_size = 0;
			m_pItems = nullptr;
		}

		Buffer( __in size_t count, __in_ecount( count ) T *pSrc ) {
			Allocate( count );
			memcpy( m_pItems, pSrc, count * sizeof( T ) );
		}

		Buffer( __in Buffer<T> &copyFrom ) {
			Allocate( copyFrom.Length() );
			memcpy( copyFrom.Ptr(), m_pItems, m_size * sizeof( T ) );
		}
	
		~Buffer () {
			if ( m_pItems != nullptr ) {
				Allocator::Release ( m_pItems );
			}
		}

		inline T* Ptr () { return m_pItems; }
	
		inline size_t Length () { return m_size; }
	
		HRESULT Allocate ( size_t length ) {
			T *temp = m_pItems;
			if ( temp != nullptr )
			{
				Allocator::Release ( temp );
			}
	
			m_pItems = Allocator::Allocate ( length );
			if ( m_pItems == nullptr )
			{
				return E_OUTOFMEMORY;
			}
			m_size = length;
	
			return S_OK;
		}
	
		T& operator[] ( size_t pos ) {
			return *( m_pItems + pos );
		}
	
	};

} }