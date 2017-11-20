using namespace RSLib;

// Setup:

void rsJustification::setJustifiedLeft()
{
  flags.setFlagTrue(justifiedLeft);
  flags.setFlagFalse(justifiedRight);
}

void rsJustification::setHorizontallyCentered()
{
  flags.setFlagFalse(justifiedLeft);
  flags.setFlagFalse(justifiedRight);
}

void rsJustification::setJustifiedRight()
{
  flags.setFlagFalse(justifiedLeft);
  flags.setFlagTrue(justifiedRight);
}

void rsJustification::setJustifiedTop()
{
  flags.setFlagTrue(justifiedTop);
  flags.setFlagFalse(justifiedBottom);
}

void rsJustification::setVerticallyCentered()
{
  flags.setFlagFalse(justifiedTop);
  flags.setFlagFalse(justifiedBottom);
}

void rsJustification::setJustifiedBottom()
{
  flags.setFlagFalse(justifiedTop);
  flags.setFlagTrue(justifiedBottom);
}

void rsJustification::setCentered()
{
  setHorizontallyCentered();
  setVerticallyCentered();
}

// Inquiry:

bool rsJustification::isJustifiedLeft() const
{
  return flags.isFlagTrue(justifiedLeft);
}

bool rsJustification::isHorizontallyCentered() const
{ 
  return flags.isFlagFalse(justifiedLeft) && flags.isFlagFalse(justifiedRight);
}

bool rsJustification::isJustifiedRight() const
{
  return flags.isFlagTrue(justifiedRight);
}

bool rsJustification::isJustifiedTop() const
{
  return flags.isFlagTrue(justifiedTop);
}

bool rsJustification::isVerticallyCentered() const
{
  return flags.isFlagFalse(justifiedTop) && flags.isFlagFalse(justifiedBottom);
}

bool rsJustification::isJustifiedBottom() const
{
  return flags.isFlagTrue(justifiedBottom);
}

bool rsJustification::isCentered() const
{
  return isHorizontallyCentered() && isVerticallyCentered();
}
