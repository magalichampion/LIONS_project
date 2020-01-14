package coolname;

import com.analog.lyric.dimple.model.domains.DiscreteDomain;
import com.analog.lyric.dimple.model.variables.Discrete;

public class Trit extends Discrete {
	
	final static DiscreteDomain domain=DiscreteDomain.range(-1,1);

	public Trit()
	{
		super(domain);
	}

}
