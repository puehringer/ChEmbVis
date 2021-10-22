import * as React from "react";
import { OverlayTrigger, Popover } from "react-bootstrap";
import { Layout } from "react-grid-layout";

export const GridItemOptions = React.memo(
  ({
    gridOptions,
    key,
    children,
    onClose,
    onSettings,
    renderInfo,
    enableMove = true,
  }: {
    gridOptions?: Partial<Layout>;
    key: string;
    children: React.ReactNode;
    onClose?(): void;
    onSettings?(): void;
    renderInfo?(): React.ReactNode;
    enableMove?: boolean;
  }) => (
    <>
      {enableMove ? (
        <i
          className="fas fa-fw fa-arrows-alt react-grid-item-hidden react-grid-item-drag-handle"
          title="Move"
          style={{
            position: "absolute",
            zIndex: 1,
            top: 3,
            left: 5,
            cursor: "pointer",
          }}
        />
      ) : null}
        {renderInfo ? (
        <OverlayTrigger
          trigger={["hover", "focus"]}
          placement="auto"
          overlay={
            <Popover id={key}>
              <Popover.Header as="h3">Projection information</Popover.Header>
              <Popover.Body>
                {renderInfo()}
              </Popover.Body>
            </Popover>
          }
        >
          <i
            className="fas fa-fw fa-info-circle react-grid-item-hidden"
            style={{
              position: "absolute",
              zIndex: 1,
              top: 3,
              left: 30,
              cursor: "pointer",
            }}
          />
        </OverlayTrigger>
      ) : null}
      {onSettings ? (
        <i
          className="fas fa-fw fa-cog react-grid-item-hidden"
          title="Change settings"
          onClick={onSettings}
          style={{
            position: "absolute",
            zIndex: 1,
            top: 3,
            left: 55,
            cursor: "pointer",
          }}
        />
      ) : null}
      {onClose ? (
        <i
          className="fas fa-fw fa-times react-grid-item-hidden"
          title="Close"
          onClick={onClose}
          style={{
            position: "absolute",
            zIndex: 1,
            top: 3,
            left: 80,
            cursor: "pointer",
          }}
        />
      ) : null}
      {children}
    </>
  )
);
