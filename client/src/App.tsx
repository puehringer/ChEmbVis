import * as React from "react";
import { InterpolationPage } from "./InterpolationPage";
import { Navbar, Container, Form, Button, Nav } from "react-bootstrap";
import { EActiveTabs, ICollection, IInterpolatedParticle, INearestNeighbors, IParticle, IRegistry } from "./interfaces";
import { EmbeddingPage } from "./EmbeddingPage";
import { downloadJSONFile, injectGetter } from "./utils";
import { ARRAY_DISTANCE_METRICS, isProxySymbol } from "./utils/constants";
import { getRegistry } from "./utils/api";
import { FileUploadModal } from "./components/FileUploadModal";
import { CollectionContext } from "./CollectionContext";

function App() {
  const [importFileModalShow, setImportFileModalShow] = React.useState<boolean>(false);
  const [activeTab, setActiveTab] = React.useState<EActiveTabs>(EActiveTabs.EMBEDDING);
  const [collections, _setCollections] = React.useState<ICollection[]>([]);
  const [registry, setRegistry] = React.useState<IRegistry | null>(null);
  const [interpolationStructures, setInterpolationStructures] = React.useState<string[]>([
    "NC1CC1C(=O)c1ccc2ccccc2c1",
    "O=C(CN1C(=O)CSc2ccc(S(=O)C3CC3)cc21)NCc1cccnc1",
  ]);

  React.useEffect(() => {
    getRegistry()
      .then((registry) => setRegistry(registry))
      .catch((e) => {
        console.error("Error initializing registry", e);
      });
  }, []);

  const setCollections = React.useCallback(
    (collections: ICollection[]) => {
      collections.forEach((collection) => {
        collection.data.forEach((p, i) => {
          p.index = i;
          p.collection = collection.name;
        });
        collection.data.forEach((p, i) => {
          // Remove jaccard distances
          Object.keys(p.properties)
            .filter((key) => key.toLowerCase().startsWith("jaccard_") || key.toLowerCase().startsWith("nn_diff_"))
            .forEach((key) => delete p.properties[key]);

          Object.entries(p.nearest_neighbors || {}).forEach(([nnKey, nn]) => {
            // Inject the actual particles into the nearest neighbors collection
            nn.knn_particles = nn.knn_ind.map((i) => collection.data[i]);
          });

          function injectEval(p: IParticle, object: any, key: string, func: string) {
            injectGetter(object, key, () => {
              try {
                return eval(func);
              } catch {
                return null;
              }
            });
          }

          function injectJaccard(object: any, key: string, embeddingName1: string, embeddingName2: string, nn: number) {
            const nn1 = p.nearest_neighbors?.[embeddingName1];
            const nn2 = p.nearest_neighbors?.[embeddingName2];
            if (nn1 && nn2) {
              injectGetter(object, key, () => {
                const xIndices = nn1.knn_ind.slice(0, nn);
                const yIndices = new Set(nn2.knn_ind.slice(0, nn));
                return xIndices.filter((x) => yIndices.has(x)).length / xIndices.length;
              });
            } else {
              throw Error(`Could not inject jaccard getter: ${embeddingName1} or ${embeddingName2} are undefined.`);
            }
          }

          function injectNNDiff(
            object: any,
            key: string,
            embeddingName: string,
            property: string,
            method: string,
            nn: number
          ) {
            const nn1 = p.nearest_neighbors?.[embeddingName];
            const value = p.properties[property];
            if (nn1 && typeof value === "number") {
              const values = nn1.knn_particles.slice(0, nn).map((p) => (p.properties[property] as number) || 0);
              injectGetter(object, key, () => ARRAY_DISTANCE_METRICS[method]?.(value, values));
            } else {
              throw Error(`Could not inject nn_diff getter: ${embeddingName} is undefined or value is not a number.`);
            }
          }

          function injectCluster(
            object: any,
            key: string,
            embeddingName: string,
          ) {
            const cluster = p.clusters?.[embeddingName];
            if (cluster) {
              injectGetter(object, key, () => cluster.label);
            } else {
              throw Error(`Could not inject cluster getter: ${embeddingName} is undefined.`);
            }
          }

          if (p.nearest_neighbors && !p.nearest_neighbors[isProxySymbol]) {
            p.nearest_neighbors = new Proxy(p.nearest_neighbors, {
              get: function (obj, prop) {
                if (prop === isProxySymbol) {
                  return true;
                }

                if (typeof prop === "symbol") {
                  return undefined;
                }

                // console.log("Nearest neighbor", obj, prop);
                if (prop.startsWith("property=") && !(prop in obj)) {
                  const [key, property] = prop.split("=");
                  const propertyValue = p.properties[property];
                  if (typeof propertyValue === "number") {
                    const lookup = new Map(
                      collection.data.map((particle) => [
                        particle,
                        Math.abs(propertyValue - (particle.properties[property]! as number)),
                      ])
                    );
                    const nearestNeighbors = collection.data
                      .sort((a, b) => lookup.get(a)! - lookup.get(b)!)
                      .slice(1, 51);

                    obj[prop] = {
                      distance_metric: "absolute_property_difference",
                      knn_dist: nearestNeighbors.map((n) => lookup.get(n)),
                      knn_ind: nearestNeighbors.map((n) => n.index),
                      knn_particles: nearestNeighbors,
                    } as INearestNeighbors;
                    // injectGetter(obj, prop, () => {})
                  } else if (typeof propertyValue === "string") {
                    const nearestNeighbors = collection.data.filter(
                      (particle) => p !== particle && particle.properties[property] === propertyValue
                    );

                    obj[prop] = {
                      distance_metric: "string_equality",
                      knn_dist: nearestNeighbors.map((n) => 0),
                      knn_ind: nearestNeighbors.map((n) => n.index),
                      knn_particles: nearestNeighbors,
                    } as INearestNeighbors;
                  }
                }
                return obj[prop];
              },
            });
          }

          if (!p.properties[isProxySymbol]) {
            p.properties = new Proxy(p.properties, {
              get: function (obj, prop) {
                if (prop === isProxySymbol) {
                  return true;
                }

                if (typeof prop === "symbol") {
                  return undefined;
                }

                try {
                  if (prop.startsWith("eval=") && !(prop in obj)) {
                    // Pattern: eval=<URL encoded function>
                    const [key, func] = prop.split("=");
                    injectEval(p, obj, prop, func);
                  }

                  if (prop.startsWith("jaccard=") && !(prop in obj)) {
                    // Pattern: jaccard=<emb1>=<emb2>=<nn>
                    const [key, emb1, emb2, nn] = prop.split("=");
                    injectJaccard(obj, prop, emb1, emb2, +nn);
                  }

                  if (prop.startsWith("nn_diff=") && !(prop in obj)) {
                    // Pattern: nn_diff=<emb>=<prop>=<method>=<nn>
                    const [key, emb, property, method, nn] = prop.split("=");
                    injectNNDiff(obj, prop, emb, property, method, +nn);
                  }

                  if (prop.startsWith("cluster=") && !(prop in obj)) {
                    // Pattern: cluster=<emb>
                    const [key, emb] = prop.split("=");
                    injectCluster(obj, prop, emb);
                  }
                } catch (e) {
                  console.error(`Error injecting ${prop} getter`, e);
                }

                return obj[prop];
              },
            });
          }
        });
      });

      _setCollections(collections);
    },
    [_setCollections]
  );

  const interpolationCollection = React.useMemo(
    () => collections.find((c) => c.name === "Interpolated") as ICollection<IInterpolatedParticle> | undefined,
    [collections]
  );

  const setInterpolationCollection = React.useCallback(
    (collection: ICollection<IInterpolatedParticle>) => {
      _setCollections((collections) => [...collections.filter((c) => c.name !== "Interpolated"), collection]);
    },
    [_setCollections]
  );

  return (
    <CollectionContext.Provider value={collections}>
      <div className="vh-100 d-flex flex-column">
        <Navbar collapseOnSelect expand="lg" bg="dark" variant="dark">
          <Container>
          <Navbar.Brand href="#home">ChEmbVis</Navbar.Brand>
          <Navbar.Toggle aria-controls="responsive-navbar-nav" />
          <Navbar.Collapse id="responsive-navbar-nav">
            <Nav className="me-auto">
              <Nav className="me-auto">
                <Nav.Link
                  href="#"
                  active={activeTab === EActiveTabs.EMBEDDING}
                  onSelect={() => setActiveTab(EActiveTabs.EMBEDDING)}
                >
                  Embedding
                </Nav.Link>
                <Nav.Link
                  href="#"
                  active={activeTab === EActiveTabs.INTERPOLATION}
                  onSelect={() => setActiveTab(EActiveTabs.INTERPOLATION)}
                >
                  Interpolation
                </Nav.Link>
                <Nav.Link href="/datasets/" target="_blank" active={false}>
                  Datasets
                </Nav.Link>
              </Nav>
            </Nav>
            <Form>
              <Button
                onClick={() => {
                  setImportFileModalShow(true);
                }}
              >
                Import
              </Button>
              <Button
                className="ms-2"
                disabled={collections.length === 0}
                onClick={() => {
                  downloadJSONFile(collections, "export");
                }}
              >
                Export
              </Button>
            </Form>
          </Navbar.Collapse>
  </Container>
        </Navbar>
        <FileUploadModal
          open={importFileModalShow}
          setOpen={setImportFileModalShow}
          onSave={(value) => {
            if (value) {
              try {
                setCollections(JSON.parse(value));
                setImportFileModalShow(false);
              } catch (e) {
                console.error("Error parsing imported file.");
              }
            }
          }}
        />
        <div className="container-fluid mt-2" style={{ flex: 1, overflow: "auto" }}>
          <div className="row" style={{ height: "100%", overflow: "auto", position: "relative" }}>
            {activeTab === EActiveTabs.EMBEDDING ? (
              <EmbeddingPage
                registry={registry}
                collections={collections}
                setCollections={setCollections}
                interpolationStructures={interpolationStructures}
                setInterpolationStructures={setInterpolationStructures}
                setActiveTab={setActiveTab}
              />
            ) : null}
            {activeTab === EActiveTabs.INTERPOLATION ? (
              <InterpolationPage
                collection={interpolationCollection}
                setCollection={setInterpolationCollection}
                structures={interpolationStructures}
                setStructures={(structures: string[]) => setInterpolationStructures(structures)}
              />
            ) : null}
          </div>
        </div>
      </div>
    </CollectionContext.Provider>
  );
}

export default App;
