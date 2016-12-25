#pragma once

namespace DataAnalysis { namespace Utils {

	// O = V * M
	template < class Numeric > void Vect4MulMat44( __in_ecount( 4 ) const Numeric *pV, __in_ecount( 16 ) const Numeric *pM, __out_ecount( 4 ) Numeric *pO ) {
		*pO = ( ( *pV ) * ( *pM ) ) + ( ( *( pV + 1 ) ) * ( *( pM + 4 ) ) ) + ( ( *( pV + 2 ) ) * ( *( pM + 8 ) ) ) + ( ( *( pV + 3 ) ) * ( *( pM + 12 ) ) );
		*( pO + 1 ) = ( ( *pV ) * ( *( pM + 1 ) ) ) + ( ( *( pV + 1 ) ) * ( *( pM + 5 ) ) ) + ( ( *( pV + 2 ) ) * ( *( pM + 9 ) ) ) + ( ( *( pV + 3 ) ) * ( *( pM + 13 ) ) );
		*( pO + 2 ) = ( ( *pV ) * ( *( pM + 2 ) ) ) + ( ( *( pV + 1 ) ) * ( *( pM + 6 ) ) ) + ( ( *( pV + 2 ) ) * ( *( pM + 10 ) ) ) + ( ( *( pV + 3 ) ) * ( *( pM + 14 ) ) );
		*( pO + 3 ) = ( ( *pV ) * ( *( pM + 3 ) ) ) + ( ( *( pV + 1 ) ) * ( *( pM + 7 ) ) ) + ( ( *( pV + 2 ) ) * ( *( pM + 11 ) ) ) + ( ( *( pV + 3 ) ) * ( *( pM + 15 ) ) );
	}

	// O = M * V
	template < class Numeric > void Mat44MulVect4( __in_ecount( 16 ) const Numeric *pM, __in_ecount( 4 ) const Numeric *pV, __out_ecount( 4 ) Numeric *pO ) {
		*pO = ( ( *pM ) * ( *pV ) ) + ( ( *( pM + 1 ) ) * ( *( pV + 1 ) ) ) + ( ( *( pM + 2 ) ) * ( *( pV + 2 ) ) ) + ( ( *( pM + 3 ) ) * ( *( pV + 3 ) ) );
		*( pO + 1 ) = ( ( *( pM + 4 ) ) * ( *( pV ) ) ) + ( ( *( pM + 5 ) ) * ( *( pV + 1 ) ) ) + ( ( *( pM + 6 ) ) * ( *( pV + 2 ) ) ) + ( ( *( pM + 7 ) ) * ( *( pV + 3 ) ) );
		*( pO + 2 ) = ( ( *( pM + 8 ) ) * ( *( pV ) ) ) + ( ( *( pM + 9 ) ) * ( *( pV + 1 ) ) ) + ( ( *( pM + 10 ) ) * ( *( pV + 2 ) ) ) + ( ( *( pM + 11 ) ) * ( *( pV + 3 ) ) );
		*( pO + 3 ) = ( ( *( pM + 12 ) ) * ( *( pV ) ) ) + ( ( *( pM + 13 ) ) * ( *( pV + 1 ) ) ) + ( ( *( pM + 14 ) ) * ( *( pV + 2 ) ) ) + ( ( *( pM + 15 ) ) * ( *( pV + 3 ) ) );
	}

} }