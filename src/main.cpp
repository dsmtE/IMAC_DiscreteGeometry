#include <boost/smart_ptr/shared_ptr.hpp>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <memory>


#include <cstdlib>
#include <cmath>
#include <iomanip>
#include "DGtal/arithmetic/LighterSternBrocot.h"
#include "DGtal/base/StdRebinders.h"

class A {};

int main() {
    // boost test 
    boost::shared_ptr<A> pA(new A);

	std::cout << pA.get() << std::endl;

	boost::shared_ptr<A> pB(pA);

	std::cout << pA.get() << std::endl;
	std::cout << pB.get() << std::endl;

    
    // DGtal
    using namespace DGtal;
    typedef DGtal::int64_t Integer;
    typedef DGtal::int64_t Quotient;
    typedef LighterSternBrocot<Integer, Quotient, StdMapRebinder> SB; // the type of the Stern-Brocot tree
    typedef SB::Fraction Fraction; // the type for fractions
    typedef std::back_insert_iterator< Fraction > OutputIterator;

    long double epsilon = 1e-14;
    long double number0 = 3.14159265;
    long double number = number0;
    ASSERT( number >= 0.0 );
    Fraction f;
    
    OutputIterator itback = std::back_inserter( f );
    Quotient i = 0;

    while ( true )
    {
        long double int_part = floorl( number );
        Quotient u = NumberTraits<long double>::castToInt64_t( int_part );
        *itback++ = std::make_pair( u, i++ );
        long double approx =
        ( (long double) NumberTraits<Integer>::castToDouble( f.p() ) )
        / ( (long double) NumberTraits<Integer>::castToDouble( f.q() ) );
        std::cout << "z = " << f.p() << " / " << f.q()
                << " =~ " << std::setprecision( 16 ) << approx << std::endl;
        number -= int_part;
        if ( ( (number0 - epsilon ) < approx )
            && ( approx < (number0 + epsilon ) ) ) break;
        number = 1.0 / number;
    }
    std::cout << "z = " << f.p() << " / " << f.q() << std::endl;
    
    return 0;
}
