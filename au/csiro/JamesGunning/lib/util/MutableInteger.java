package au.csiro.JamesGunning.lib.util;

/**
   Simple wrapper class of a basic integer that is mutable, and
   thus can be reused to save garbage collection burden for
   Collections etc that use a lot of Integers.

   Should also have all the methods of Integer(), but these remain
   to be implemented.
*/

public final class MutableInteger
{
  public int value;
  
  /**
	 Constructor: set the value to that passed in
  */
  public MutableInteger(int value)
	{
	  this.value = value;
	}
  
  public MutableInteger()
	{
	}

  /**
	 Chamges the value to a new value
  */
  public MutableInteger set(int value)
	{
	  this.value = value;
	  return this;
	}

  public int intValue()
  {
	  return value;
  }

  // to use this class as keys in a Hashmap, we need:
  public int hashCode()
  {
	 return value;
  }

  // to use this class as keys in a Hashmap, we need:
  public boolean equals(Object obj)
  {
	if (obj == null) return false;

	if(getClass() == obj.getClass())  
	  {
		if (value == ((MutableInteger)obj).value)
		  {
			return true;
		  }
		else
		  return false;
	  }
	else
	  return false;
  }
}
