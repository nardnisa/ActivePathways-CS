package org.jurilab.apw.internal.demo5;
import static org.cytoscape.work.ServiceProperties.IN_MENU_BAR;
import static org.cytoscape.work.ServiceProperties.PREFERRED_MENU;
import static org.cytoscape.work.ServiceProperties.TITLE;
import java.util.Properties;
import org.cytoscape.work.TaskFactory;
import org.cytoscape.service.util.CyServiceRegistrar;
import org.cytoscape.view.model.CyNetworkViewManager;
import org.cytoscape.view.vizmap.VisualMappingManager;
import org.cytoscape.service.util.AbstractCyActivator;
import org.osgi.framework.BundleContext;
import org.osgi.framework.InvalidSyntaxException;

public class CyActivator extends AbstractCyActivator {
	public CyActivator() {
		super();
	}
	public void start(BundleContext bc) throws InvalidSyntaxException {
		CyServiceRegistrar serviceRegistrar = getService(bc, CyServiceRegistrar.class);
		CyNetworkViewManager viewManager = getService(bc, CyNetworkViewManager.class);
		viewManager.getNetworkViewSet();
		APWManager manager = new APWManager(serviceRegistrar);
		VisualMappingManager vizManager = manager.getService(VisualMappingManager.class);
		 
		
		
		Properties runActivePathways = new Properties();
		CreateAPWTaskFactory createAPWTaskFactory = new CreateAPWTaskFactory(manager);
		runActivePathways.setProperty(PREFERRED_MENU, "Apps");
		runActivePathways.setProperty(TITLE, "ActivePathways");
		runActivePathways.setProperty(IN_MENU_BAR,"true");
		registerService(bc, createAPWTaskFactory, TaskFactory.class, runActivePathways);
	}
}

