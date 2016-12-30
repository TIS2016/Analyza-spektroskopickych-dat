//
//  Buffer.h
//  Utils
//
//  Created by Filip Kovac on 19/10/2016.
//  Copyright © 2016 Filip Kovac. All rights reserved.
//

#pragma once

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

		Buffer( __in const size_t count, __in_ecount( count ) const T *pSrc ) {
			Set( count, pSrc );
		}

		Buffer( __in const Buffer<T> &copyFrom ) {
			Set( copyFrom );
		}
	
		~Buffer () {
			if ( m_pItems != nullptr ) {
				Allocator::Release ( m_pItems );
			}
		}

		inline T* Ptr () const { return m_pItems; }
	
		inline size_t Length () const { return m_size; }
	
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

		HRESULT Set( __in const size_t count, __in_ecount(count) const T *pSrc ) {
			if ( pSrc == nullptr ) {
				return E_INVALIDARG;
			}
			
			HRESULT hr = Allocate( count );
			if ( SUCCEEDED( hr ) ) {
				// for safety reasons, do copy in loop
				T *pDst = m_pItems;
				for ( size_t i = 0; i < count; i++ ) {
					*pDst = *pSrc;
					pDst++;
					pSrc++;
				}
			}
			
			return hr;
		}

		HRESULT Set( __in const Buffer<T> &src ) {
			return Set( src.Length(), src.Ptr() );
		}

		HRESULT Set( __in const vector<T> &src ) {
			return Set( src.size(), src.data() );
		}
	
		T& operator[] ( size_t pos ) const {
			return *( m_pItems + pos );
		}
	
	};

} }