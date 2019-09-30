package org.jurilab.apw.internal.demo5;

import org.cytoscape.work.AbstractTaskFactory;
import org.cytoscape.work.TaskIterator;

public class CreateAPWTaskFactory extends AbstractTaskFactory {
	public final APWManager manager;

	public CreateAPWTaskFactory(final APWManager manager) {
		this.manager = manager;
	}
	@Override
	public TaskIterator createTaskIterator() {
		TaskIterator ti = new TaskIterator(new CreateAPWTask(manager));
		return ti;
	}
	
}
