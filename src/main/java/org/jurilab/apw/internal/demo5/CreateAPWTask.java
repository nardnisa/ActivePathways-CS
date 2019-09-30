package org.jurilab.apw.internal.demo5;
import org.cytoscape.work.AbstractTask;
import org.cytoscape.work.TaskMonitor;

public class CreateAPWTask extends AbstractTask {
	final APWManager manager;
	public CreateAPWTask(final APWManager manager) {
		this.manager = manager;
	}
	public void run(TaskMonitor taskMonitor) {
		taskMonitor.setTitle("Creating Enrichment Map");
		manager.initComponents_simple_initial();
	}
}

