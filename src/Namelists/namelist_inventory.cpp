#include "copyright.h"
#include "namelist_inventory.h"

namespace stormm {
namespace namelist {

//-------------------------------------------------------------------------------------------------
NamelistToken::NamelistToken(const std::string &title_in, NmlFuncPtr producer_in) :
    title{title_in}, producer{producer_in}
{}

//-------------------------------------------------------------------------------------------------
const std::string& NamelistToken::getTitle() const {
  return title;
}

//-------------------------------------------------------------------------------------------------
NamelistEmulator NamelistToken::invoke(const TextFile &tf, int *start_line, bool *found,
                                       const ExceptionResponse policy,
                                       const WrapTextSearch wrap) const {
  return producer(tf, start_line, found, policy, wrap);
}

} // namespace namelist
} // namespace stormm
